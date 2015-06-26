#################################################################################################
## Fitting procedure for multivariate mixtures of Erlangs (MME) to censored and truncated data ##
#################################################################################################

## Initial values

# supply the maximum initial number of shapes M in each dimension and the spread factor

.MME_initial <- function(lower, upper, trunclower=rep(0, ncol(lower)), truncupper=rep(Inf, ncol(lower)), M=10, s=100){
  n <- nrow(lower)
  d <- ncol(lower)
  # initializing data (uncensored + lower bound for right censored data points)
  initial_data <- lower
  # take upper bound for left censored data points
  initial_data[is.na(lower)] <- upper[is.na(lower)]
  # take mean for interval censored data points
  initial_data[lower!=upper & !is.na(lower!=upper)] <- (lower[lower!=upper & !is.na(lower!=upper)] + upper[lower!=upper & !is.na(lower!=upper)]) / 2
  # replace 0 initial data (= right censored at 0) by NA, since they don't add any information and will cause initial shape = 0
  initial_data[initial_data == 0] <- NA
  # non-observed values remain NA in initial_data
  # initialize shapes
  maxima <- apply(initial_data, 2, max, na.rm = TRUE)
  theta <-  min(maxima) / s
  r <- vector("list", d)
  for(j in 1:d){
    # use quantiles for initial shape choices
    r[[j]] <- unique(ceiling(quantile(initial_data[, j], probs = seq(0, 1, length.out = M), na.rm = TRUE)/theta))
  }
  # all shape combinations
  shape <- as.matrix(expand.grid(r))
  # initialize weights
  alpha <- rep(0, nrow(shape))
  indicator <- vector(mode = "list", length = length(r))
  for(i in 1:n){
    for(j in 1:length(r)){
      indicator[[j]][1] <- (initial_data[i, j] <= r[[j]][1]*theta)
      for(k in 2:length(r[[j]])){
        indicator[[j]][k] <- ( r[[j]][k-1]*theta < initial_data[i, j] & initial_data[i, j]  <= r[[j]][k]*theta)
      }
      indicator[[j]][is.na(indicator[[j]])] <- 1/length(r[[j]])
    }
    alpha <- alpha + rowProds(as.matrix(expand.grid(indicator)))
  }
  shape <- shape[alpha>0, , drop = FALSE]
  alpha <- alpha[alpha>0]/sum(alpha)
  # alpha to beta
  t_probabilities <- colProds(matrix( pgamma(truncupper, shape=t(shape), scale=theta) - pgamma(trunclower, shape=t(shape), scale=theta), dim(shape)[2], dim(shape)[1]))
  beta <- alpha * t_probabilities / sum(alpha*t_probabilities)
  list(theta=theta, shape=shape, alpha=alpha, beta=beta)
}

## Log likelihood

.MME_loglikelihood <- function(f, beta, t_probabilities){
  likelihood_contribution <- rowSums(t(t(f)*beta/t_probabilities))
  loglikelihood_contribution <- ifelse(likelihood_contribution>0, log(likelihood_contribution), -1000)
  # loglikelihood
  sum(loglikelihood_contribution)
}

## multivariate Erlang density/probability for each observation and shape combination

.MME_f <- function(lower, upper, shape, theta){
  f = matrix(0, nrow = nrow(lower), ncol = nrow(shape))
  for(j in 1:nrow(shape)){
    f[,j] <- rowProds(ifelse(lower == upper, t(dgamma(t(lower), shape=shape[j,], scale=theta)), t(pgamma(t(upper), shape=shape[j,], scale=theta)) - t(pgamma(t(lower), shape=shape[j,], scale=theta))))
  }
  f
}

## z_{ir}^{(k)}: posterior probabilities

.Estep_z <- function(f, beta, t_probabilities, R){
  components <- t(t(f)*beta/t_probabilities)
  z <- components / rowSums(components)
  # in case all z_{ir} for all r are numerically 0
  z[is.nan(z)] <- 1/R
  z
}

## Expected value of the sum of the elements of an observation

.Estep_EX <-function(lower, upper, shape, theta, z){
  EX_r <- matrix(0, nrow = nrow(lower), ncol = nrow(shape))
  for(j in 1:nrow(shape)){
    EX_r_ind <- ifelse(lower == upper, lower, t(shape[j,]*theta*(pgamma(t(upper), shape=shape[j,]+1, scale=theta) - pgamma(t(lower), shape=shape[j,]+1, scale=theta))/(pgamma(t(upper), shape=shape[j,], scale=theta) - pgamma(t(lower), shape=shape[j,], scale=theta))))
    # replace numerical 0/0 (NaN) or Inf by correct expected value
    EX_r_ind <- ifelse(is.nan(EX_r_ind) | EX_r_ind==Inf, ifelse(t(t(lower) > shape[j,]*theta), lower, upper), EX_r_ind)
    EX_r[,j] <- rowSums(EX_r_ind)
  }
  EX <- rowSums(EX_r*z)
}

## Correction term for specific dimension (auxiliary function)

.Mstep_T_j <- function(trunclower, truncupper, shape, theta){
  # avoid NaN
  if(truncupper==Inf){
    # take log first for numerical stability (avoid Inf / Inf)
    deriv_trunc_log_1 <- shape*log(trunclower)-trunclower/theta - (shape-1)*log(theta) - lgamma(shape) - log(1 - pgamma(trunclower, shape, scale=theta))
    deriv_trunc <- exp(deriv_trunc_log_1)
  } else{
    deriv_trunc_log_1 <- shape*log(trunclower)-trunclower/theta - (shape-1)*log(theta) - lgamma(shape) - log(pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta))
    deriv_trunc_log_2 <- shape*log(truncupper)-truncupper/theta - (shape-1)*log(theta) - lgamma(shape) - log(pgamma(truncupper, shape, scale=theta) - pgamma(trunclower, shape, scale=theta))
    deriv_trunc <- exp(deriv_trunc_log_1)-exp(deriv_trunc_log_2)
  }
  deriv_trunc
}

## T^{(k)}: correction term due to truncation

.Mstep_T <- function(trunclower, truncupper, shape, theta, beta){
  t <- matrix(0, nrow = nrow(shape), ncol = ncol(shape))
  for(j in 1:ncol(shape)){
      t[,j] = .Mstep_T_j(trunclower[j], truncupper[j], shape[,j], theta)
  }
  sum(beta*rowSums(t))
}

## Auxiliary function used to maximize theta in the M-step

.theta_nlm <- function(theta, EX, n, beta, shape, trunclower, truncupper){
  T <- .Mstep_T(trunclower, truncupper, shape, theta, beta)
  (theta - (sum(EX)/n-T)/sum(beta*rowSums(shape)))^2
}

## EM algorithm

.MME_em <- function(lower, upper, trunclower=rep(0, ncol(lower)), truncupper=rep(Inf, ncol(lower)), theta, shape, beta, eps=1e-03, print=TRUE, beta_tol = 10^(-5), max_iter = Inf){
  lower <- as.matrix(lower)
  upper <- as.matrix(upper)
  n <- nrow(lower)
  d <- ncol(lower)
  R <- nrow(shape)
  for(j in 1:d){
    lower[,j][is.na(lower[,j])] = trunclower[j]
    upper[,j][is.na(upper[,j])] = truncupper[j]
  }
  iteration <- 1
  # multivariate Erlang density/probability for each observation and shape combination
  f <- .MME_f(lower, upper, shape, theta)
  # multivariate Erlang truncations probabilities
  t_probabilities <- colProds(matrix( pgamma(truncupper, shape=t(shape), scale=theta) - pgamma(trunclower, shape=t(shape), scale=theta) , dim(shape)[2], dim(shape)[1]))
  loglikelihood <- .MME_loglikelihood(f, beta, t_probabilities)
  old_loglikelihood <- -Inf
  history_loglikelihood <- loglikelihood
  while(loglikelihood - old_loglikelihood > eps & iteration < max_iter){
    old_loglikelihood <- loglikelihood
    # E step
    z <- .Estep_z(f, beta, t_probabilities, R)
    EX <- .Estep_EX(lower, upper, shape, theta, z)
    # M step
    beta <- colSums(z)/n
    if(min(beta) < beta_tol){
      shape <- shape[beta > beta_tol, , drop=FALSE]
      R <- nrow(shape)
      beta <- beta[beta > beta_tol]
      beta <- beta/sum(beta)
    }
    theta <- nlm(.theta_nlm, theta, EX, n, beta, shape, trunclower, truncupper)$estimate
    iteration <- iteration + 1
    # multivariate Erlang density/probability for each observation and shape combination
    f <- .MME_f(lower, upper, shape, theta)
    # multivariate Erlang truncations probabilities
    t_probabilities <- colProds(matrix( pgamma(truncupper, shape=t(shape), scale=theta) - pgamma(trunclower, shape=t(shape), scale=theta) , dim(shape)[2], dim(shape)[1]))
    loglikelihood <- .MME_loglikelihood(f, beta, t_probabilities)
    if(print) print(loglikelihood)
    history_loglikelihood <- c(history_loglikelihood, loglikelihood)
  }
  # beta to alpha
  alpha_tilde <- beta / t_probabilities
  alpha <- alpha_tilde / sum(alpha_tilde)
  list(alpha = alpha, beta = beta, shape = shape, theta = theta, R=R, loglikelihood = loglikelihood, history_loglikelihood = history_loglikelihood, iteration = iteration, AIC=-2*loglikelihood+2*(length(alpha)*(ncol(lower)+1)),BIC=-2*loglikelihood+log(nrow(lower))*(length(alpha)*(ncol(lower)+1)))
}

## Shape adjustments

.MME_shape_adj <- function(lower, upper, trunclower, truncupper, theta, shape, beta, eps=1e-03, print=TRUE, beta_tol = 10^(-5), max_iter = Inf){
  fit <- .MME_em(lower, upper, trunclower, truncupper, theta, shape, beta, eps, print=FALSE, beta_tol, max_iter)
  loglikelihood <- fit$loglikelihood
  shape <- fit$shape
  theta <- fit$theta
  beta <- fit$beta
  alpha <- fit$alpha
  R <- dim(shape)[1]
  d <- dim(shape)[2]
  # before and after are the loglikelihoods used in the outer while loop
  before_loglikelihood <- -Inf
  after_loglikelihood <- loglikelihood
  while(after_loglikelihood > before_loglikelihood + eps){
    before_loglikelihood <- after_loglikelihood
    # Loop over dimensions
    for(j in 1:d){
      # Try increasing the shapes
      i <- R
      while(i>0){
        repeat{
          new_shape <- shape
          new_shape[i,j] <- new_shape[i,j]+1
          # Test for new shape combination
          if( any(colAlls(t(shape) == new_shape[i,])) ) break
          fit <- .MME_em(lower, upper, trunclower, truncupper, theta, new_shape, beta, eps, print=FALSE, beta_tol, max_iter)
          new_loglikelihood <- fit$loglikelihood
          if(new_loglikelihood > loglikelihood + eps){
            loglikelihood <- new_loglikelihood
            shape <- fit$shape
            theta <- fit$theta
            beta <- fit$beta
            alpha <- fit$alpha
            # number of shape combinations might have changed after EM algorithm if beta_tol > 0
            R <- dim(shape)[1]
            i <- min(i, R)
            if(print) {
              cat("loglikelihood = ", loglikelihood, "theta = ", theta, "\n", "alpha, shape = ", "\n")
              print(cbind(alpha = alpha, shape))
            }
          } else break
        }
        i <- i-1
      }
      # Try decreasing the shapes
      i <- 1
      while(i <= R){
        repeat{
          new_shape <- shape
          new_shape[i,j] <- new_shape[i,j]-1
          if( any(colAlls(t(shape) == new_shape[i,])) | new_shape[i,j]==0 ) break
          fit <- .MME_em(lower, upper, trunclower, truncupper, theta, new_shape, beta, eps, print=FALSE, beta_tol, max_iter)
          new_loglikelihood <- fit$loglikelihood
          if(new_loglikelihood > loglikelihood + eps){
            loglikelihood <- new_loglikelihood
            shape <- fit$shape
            theta <- fit$theta
            beta <- fit$beta
            alpha <- fit$alpha
            # number of shape combinations might have changed after EM algorithm if beta_tol > 0
            R <- dim(shape)[1]
            i <- min(i, R)
            if(print) {
              cat("loglikelihood = ", loglikelihood, "theta = ", theta, "\n", "alpha, shape = ", "\n")
              print(cbind(alpha = alpha, shape))
            }
          } else break
        }
        i <- i + 1
      }
    }
    after_loglikelihood <- loglikelihood
  }
  list(alpha = alpha, beta = beta, shape = shape, theta = theta, R=R, loglikelihood = loglikelihood, AIC=-2*loglikelihood+2*(length(alpha)*(ncol(lower)+1)),BIC=-2*loglikelihood+log(nrow(lower))*(length(alpha)*(ncol(lower)+1)))
}

## Shape reduction

.MME_shape_red <- function(lower, upper, trunclower=rep(0, ncol(lower)), truncupper=rep(Inf, ncol(lower)), theta, shape, beta, criterium="BIC", eps=1e-03, print=TRUE, beta_tol = 10^(-5), max_iter = Inf, adj = TRUE){
  if( adj == TRUE ){
    fit <- .MME_shape_adj(lower, upper, trunclower, truncupper, theta, shape, beta, eps, print, beta_tol, max_iter)
  } else{
    fit <- .MME_em(lower, upper, trunclower, truncupper, theta, shape, beta, eps, print, beta_tol, max_iter)
  }
  loglikelihood <- fit$loglikelihood
  IC <- fit[[criterium]]
  shape <- fit$shape
  theta <- fit$theta
  beta <- fit$beta
  alpha <- fit$alpha
  R <- dim(shape)[1]
  if(print) cat("Number of shape combinations = ", R, ", ", criterium, " = ", IC, "\n")
  improve <- TRUE
  while((improve==TRUE) && R > 1){
    new_shape <- shape[beta != min(beta), , drop=FALSE]
    new_beta <- beta[beta != min(beta)]
    new_beta <- new_beta/sum(new_beta)
    if( adj == TRUE){
      fit <- .MME_shape_adj(lower, upper, trunclower, truncupper, theta, new_shape, new_beta, eps, print, beta_tol, max_iter)
    } else{
      fit <- .MME_em(lower, upper, trunclower, truncupper, theta, new_shape, new_beta, eps, print, beta_tol, max_iter)
    }
    new_IC <- fit[[criterium]]
    if(new_IC < IC){
      IC <- new_IC
      loglikelihood <- fit$loglikelihood
      shape <- fit$shape
      theta <- fit$theta
      beta <- fit$beta
      alpha <- fit$alpha
      R <- dim(shape)[1]
      if(print) cat("Number of shape combinations = ", R, ", ", criterium, " = ", IC, "\n")
    } else{improve <- FALSE}
  }
  list(alpha = alpha, beta = beta, shape = shape, theta = theta, R=R, loglikelihood = loglikelihood, AIC=-2*loglikelihood+2*(length(alpha)*(ncol(lower)+1)),BIC=-2*loglikelihood+log(nrow(lower))*(length(alpha)*(ncol(lower)+1)))
}

## Calibration procedure for multivariate mixtures of Erlangs by repeatedly using the EM algorithm while reducing the shape parameters based on an information criterium (AIC and BIC implemented)
## First reduction of the shape combinations, afterwards adjustment and reduction of the shape combinations.

.MME_fit <- function(lower, upper = lower, trunclower = rep(0, ncol(as.matrix(lower))), truncupper = rep(Inf, ncol(as.matrix(lower))), M=10, s=100, criterium="BIC", eps=1e-03, print=TRUE, beta_tol = 10^(-5), max_iter = Inf, initial){
  lower <- as.matrix(lower)
  upper <- as.matrix(upper)
  if(missing(initial)) initial <- .MME_initial(lower, upper, trunclower, truncupper, M, s)
  # Reduction of the shape combinations
  fit_red <- .MME_shape_red(lower, upper, trunclower, truncupper, initial$theta, initial$shape, initial$beta, criterium, eps, print, beta_tol, max_iter, adj = FALSE)
  # Subsequent adjustment and reduction of the shape combinations
  fit_adj <- .MME_shape_red(lower, upper, trunclower, truncupper, fit_red$theta, fit_red$shape, fit_red$beta, criterium, eps, print, beta_tol, max_iter, adj = TRUE)
  list(alpha = fit_adj$alpha, beta = fit_adj$beta, shape = fit_adj$shape, theta = fit_adj$theta, R=fit_adj$R, loglikelihood = fit_adj$loglikelihood, AIC=fit_adj$AIC, BIC=fit_adj$BIC, M=M, s=s, alpha_red = fit_red$alpha, beta_red = fit_red$beta, shape_red = fit_red$shape, theta_red = fit_red$theta, R_red=fit_red$R, loglikelihood_red = fit_red$loglikelihood, AIC_red=fit_red$AIC, BIC_red=fit_red$BIC, alpha_initial = initial$alpha, beta_initial = initial$beta, theta_initial = initial$theta, shape_initial = initial$shape)
}



#' @title Fits univariate and multivariate mixtures of Erlang distributions to possibly censored and/or truncated data.
#'
#' @param lower A matrix specifying the lower censoring points, observations in rows, dimensions in coloms.
#' @param upper A matrix specifying the upper censoring points, observations in rows, dimensions in coloms.
#' @param trunclower A vector specifying the lower truncation points in each dimension.
#' @param truncupper A vector specifying the upper truncation points in each dimension.
#' @param M A vector of values for the tuning parameter M.
#' @param s A vector of values for the tuning parameter s (the spread).
#' @param nCores Number of cores available for parallel computation.
#' @param criterium Character vector specifying information criterium to use, either "AIC" or "BIC".
#' @param eps Covergence threshold used in the EM algorithm.
#' @param beta_tol Threshold for the mixing weights below which the corresponding shape parameter vector is considered neglectable.
#' @param print Bool: print intermediate results, either TRUE or FALSE.
#' @param file If print is TRUE, specify file name.
#' @param max_iter Maximum number of iterations in a single EM algorithm.
#'
#' @return A list containing:
#' best_model: The final MME, judged to be the best according to the criterium used. Value is a list.
#' performances: A matrix summarizing the performance of the different fitted MME models (M, s, criterium, R).
#' all_model: A list containing all fitted MME models.
#' @examples
#' \dontrun{
#' library(MASS)
#' MME_tune(lower = geyser, upper = geyser, trunclower = c(0, 0),
#' truncupper = c(Inf, Inf), M = c(5, 10, 20), s = c(seq(10, 100, 10), 200),
#' nCores = detectCores(), criterium = "BIC", eps = 1e-03,
#'  beta_tol = 10^(-5), print=TRUE, file="log.txt", max_iter = Inf)
#' }
MME_tune <- function(lower, upper = lower, trunclower = rep(0, ncol(as.matrix(lower))), truncupper = rep(Inf, ncol(as.matrix(lower))), M = 10, s = 100, nCores = detectCores(), criterium = "BIC", eps = 1e-03, beta_tol = 10^(-5), print=TRUE, file="log.txt", max_iter = Inf){
  tuning_parameters <- expand.grid(M, s)
  cl <- makePSOCKcluster(nCores)
  clusterExport(cl, c("lower", "upper", "trunclower", "truncupper", "tuning_parameters", "criterium", "eps", "print", "beta_tol", "max_iter", ".MME_initial", ".MME_loglikelihood", ".MME_f", ".Estep_z", ".Estep_EX", ".Mstep_T_j", ".Mstep_T", ".theta_nlm", ".MME_em", ".MME_shape_adj", ".MME_shape_red", ".MME_fit"), env = environment())
  registerDoParallel(cl)
  if(print) writeLines(c(""), file)
  all_model <- foreach(i = 1:nrow(tuning_parameters), .packages='matrixStats', .errorhandling = 'remove') %dopar% {
    if(print) cat(paste("M = ", tuning_parameters[i, 1], ", s = ", tuning_parameters[i, 2], "\n"), file = file, append = TRUE)
    .MME_fit(lower, upper, trunclower, truncupper, M = tuning_parameters[i, 1], s = tuning_parameters[i, 2], criterium, eps, FALSE, beta_tol, max_iter)
  }
  stopCluster(cl)
  performances <- data.frame(sapply(all_model, with, M), sapply(all_model, with, s),  sapply(all_model, function(x) with(x, get(criterium))), sapply(all_model, with, R))
  colnames(performances) <- c('M', 's', criterium, 'R')
  best_index <- which.min(performances[, criterium])
  best_model <- all_model[[best_index]]
  list(best_model = best_model, performances = performances, all_model = all_model)
}
