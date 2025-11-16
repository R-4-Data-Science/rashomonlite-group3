fit_models_dim <- function(X, y, k, m) {
  all_vars <- colnames(X)
  n_vars <- length(all_vars)

  results <- vector("list", m)
  n_generated <- 0
  seen <- character(0)
  max_iterations <- m * 20
  iterations <- 0

  #k = 1 models
  if (k == 1) {
    n_to_fit <- min(m, n_vars)
    vars_to_use <- if (n_vars > m) sample(all_vars, m) else all_vars

    for (i in seq_along(vars_to_use)) {
      vars_i <- vars_to_use[i]
      formula_i <- as.formula(paste("y ~", vars_i))
      model_i <- glm(formula_i, data = data.frame(y = y, X), family = binomial)

      coef_i <- coef(model_i)[-1]
      pvals_i <- summary(model_i)$coefficients[-1, 4]

      results[[i]] <- list(
        vars = vars_i,
        aic = AIC(model_i),
        coef = coef_i,
        pvals = pvals_i
      )
      n_generated <- n_generated + 1
    }
  } else {
    #k > 1 models
    while (n_generated < m && iterations < max_iterations) {
      iterations <- iterations + 1

      #randomly sample k predictors
      vars_i <- sample(all_vars, k)

      #create signature to check uniqueness
      sig <- paste(sort(vars_i), collapse = ",")

      if (sig %in% seen) next

      #mark as seen
      seen <- c(seen, sig)
      n_generated <- n_generated + 1

      #fit model
      formula_i <- as.formula(paste("y ~", paste(vars_i, collapse = " + ")))
      model_i <- glm(formula_i, data = data.frame(y = y, X), family = binomial)

      #extract coeff and pvals
      coef_i <- coef(model_i)[-1]
      pvals_i <- summary(model_i)$coefficients[-1, 4]

      #store result
      results[[n_generated]] <- list(
        vars = vars_i,
        aic = AIC(model_i),
        coef = coef_i,
        pvals = pvals_i
      )
    }
  }

  #trim to actual number generated
  results <- results[1:n_generated]

  #convert to data frame
  df <- data.frame(
    vars = I(lapply(results, function(x) x$vars)),
    aic = sapply(results, function(x) x$aic),
    coef = I(lapply(results, function(x) x$coef)),
    pvals = I(lapply(results, function(x) x$pvals)),
    stringsAsFactors = FALSE
  )

  return(df)
}

#------------------------------------
#select top alpha% of models by AIC
select_best <- function(models_df, alpha = 0.5) {
  sorted_df <- models_df[order(models_df$aic), ]
  n_keep <- ceiling(nrow(sorted_df) * alpha)
  return(sorted_df[1:n_keep, ])
}

#------------------------------------

grow_from <- function(best_df, all_vars, k_next, m, X, y) {
  results <- vector("list", m)
  n_generated <- 0
  seen <- character(0)

  max_iterations <- m * 20
  iterations <- 0

  while (n_generated < m && iterations < max_iterations) {
    iterations <- iterations + 1

    #sample a base model from best_df
    base_idx <- sample(nrow(best_df), 1)
    base_vars <- best_df$vars[[base_idx]]

    #find variables not in  model
    available_vars <- setdiff(all_vars, base_vars)
    if (length(available_vars) == 0) next

    #sample one new variable
    new_var <- sample(available_vars, 1)
    new_vars <- sort(c(base_vars, new_var))

    #create signature to check uniqueness
    sig <- paste(new_vars, collapse = ",")
    if (sig %in% seen) next  # Skip duplicates

    #mark as seen
    seen <- c(seen, sig)
    n_generated <- n_generated + 1

    #fit the model
    formula_i <- as.formula(paste("y ~", paste(new_vars, collapse = " + ")))
    model_i <- glm(formula_i, data = data.frame(y = y, X), family = binomial)

    #extract coeff and pvals
    coef_i <- coef(model_i)[-1]
    pvals_i <- summary(model_i)$coefficients[-1, 4]

    #store result
    results[[n_generated]] <- list(
      vars = new_vars,
      aic = AIC(model_i),
      coef = coef_i,
      pvals = pvals_i
    )
  }

  #trim
  results <- results[1:n_generated]

  #convert to data frame
  df <- data.frame(
    vars = I(lapply(results, function(x) x$vars)),
    aic = sapply(results, function(x) x$aic),
    coef = I(lapply(results, function(x) x$coef)),
    pvals = I(lapply(results, function(x) x$pvals)),
    stringsAsFactors = FALSE
  )

  return(df)
}

#-----------------------------------

#' Run Rashomon Model Expansion Procedure
#'
#' Executes a sequential model–growing workflow to explore a Rashomon set of
#' models across increasing model dimensions. For each dimension \eqn{k}, the
#' function fits many candidate models, evaluates them using an AIC-based
#' selection rule, and retains only the near-optimal subset. These retained
#' models are expanded to dimension \eqn{k+1}, continuing until \code{pmax}.
#' The function also produces summaries of predictor usage and AIC distributions.
#'
#' @param X A data frame or matrix containing predictor variables. Column names
#'   are treated as variable identifiers.
#' @param y A numeric response vector with length equal to the number of rows in
#'   \code{X}.
#' @param pmax Integer. Maximum model dimension (largest number of predictors
#'   allowed in any model). Defaults to \code{5}.
#' @param m Integer. Number of candidate models to generate at each dimension.
#'   Passed to internal functions such as \code{fit_models_dim()} and
#'   \code{grow_from()}. Defaults to \code{90}.
#' @param alpha Numeric. AIC tolerance parameter used by \code{select_best()} to
#'   identify near-optimal models. Smaller values are more selective. Defaults
#'   to \code{0.5}.
#'
#' @details
#' The procedure follows these steps:
#' \enumerate{
#'   \item Initializes storage for models across dimensions.
#'   \item Fits all 1-predictor models using \code{fit_models_dim()} and selects
#'         the best subset with \code{select_best()}.
#'   \item For dimensions \eqn{k = 2, \ldots, pmax}, expands the previously
#'         selected models using \code{grow_from()} and again filters them with
#'         \code{select_best()}.
#'   \item Computes counts of how often each predictor appears in the selected
#'         models at each dimension.
#'   \item Produces an AIC summary table useful for boxplots and other
#'         diagnostics.
#' }
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{M}}{A list of length \code{pmax}. Each element \code{M[[k]]}
#'     contains a data frame of the best retained models at dimension \eqn{k},
#'     typically including fields such as \code{vars} (predictor sets) and
#'     \code{aic}.}
#'
#'   \item{\code{predictor_counts}}{A list of data frames. For each dimension,
#'     reports how many times each predictor appears among the retained models.}
#'
#'   \item{\code{aic_summary}}{A data frame with columns \code{dimension} and
#'     \code{aic}, containing the AIC values for all retained models across all
#'     dimensions.}
#' }
#'
#' @seealso
#' \code{\link{fit_models_dim}}, \code{\link{grow_from}}, \code{\link{select_best}}
#' @author Jonah Kennedy
#' @export
#' @examples
#' result <- run_rashomon(X, y, pmax = 4, m = 100, alpha = 0.3)
#' result$aic_summary
#' result$predictor_counts$k2
run_rashomon <- function(X, y, pmax = 5, m = 90, alpha = 0.5) {
  all_vars <- colnames(X)

  #Initialization
  M <- vector("list", pmax)
  names(M) <- paste0("k", 1:pmax)

  #fit all 1st dim models
  cat("Fitting models for k = 1...\n")
  models_k1 <- fit_models_dim(X, y, k = 1, m = m)
  best_k1 <- select_best(models_k1, alpha = alpha)
  M[[1]] <- best_k1

  #fit 2-p_max models
  for (k in 2:pmax) {
    cat("Fitting models for k =", k, "...\n")
    models_k <- grow_from(M[[k-1]], all_vars, k_next = k, m = m, X = X, y = y)
    best_k <- select_best(models_k, alpha = alpha)
    M[[k]] <- best_k
  }

  #count of models per predictor at each dimension
  predictor_counts <- vector("list", pmax)
  names(predictor_counts) <- paste0("k", 1:pmax)

  for (k in 1:pmax) {
    counts <- table(unlist(M[[k]]$vars))
    predictor_counts[[k]] <- data.frame(
      predictor = names(counts),
      count = as.numeric(counts),
      stringsAsFactors = FALSE
    )
  }

  #AIC summary data for boxplots
  aic_summary <- data.frame(
    dimension = rep(1:pmax, sapply(M, nrow)),
    aic = unlist(lapply(M, function(df) df$aic)),
    stringsAsFactors = FALSE
  )

  #return results
  return(list(
    M = M,
    predictor_counts = predictor_counts,
    aic_summary = aic_summary
  ))
}
