#################################################################################################
## Fitting procedure for multivariate mixtures of Erlangs (MME) to censored and truncated data ##
#################################################################################################

#' @import foreach 
NULL

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
    alpha <- alpha + matrixStats::rowProds(as.matrix(expand.grid(indicator)))
  }
  shape <- shape[alpha>0, , drop = FALSE]
  alpha <- alpha[alpha>0]/sum(alpha)
  # alpha to beta
  t_probabilities <- matrixStats::colProds(matrix( pgamma(truncupper, shape=t(shape), scale=theta) - pgamma(trunclower, shape=t(shape), scale=theta), dim(shape)[2], dim(shape)[1]))
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
    f[,j] <- matrixStats::rowProds(ifelse(lower == upper, t(dgamma(t(lower), shape=shape[j,], scale=theta)), t(pgamma(t(upper), shape=shape[j,], scale=theta)) - t(pgamma(t(lower), shape=shape[j,], scale=theta))))
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
  t_probabilities <- matrixStats::colProds(matrix( pgamma(truncupper, shape=t(shape), scale=theta) - pgamma(trunclower, shape=t(shape), scale=theta) , dim(shape)[2], dim(shape)[1]))
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
    t_probabilities <- matrixStats::colProds(matrix( pgamma(truncupper, shape=t(shape), scale=theta) - pgamma(trunclower, shape=t(shape), scale=theta) , dim(shape)[2], dim(shape)[1]))
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
          if( any(matrixStats::colAlls(t(shape) == new_shape[i,])) ) break
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
          if( any(matrixStats::colAlls(t(shape) == new_shape[i,])) | new_shape[i,j]==0 ) break
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

# Check input arguments for MEtune
.ME_checkInput  <- function(lower, upper, trunclower, truncupper) {
  
  nl <- length(lower)
  nu <- length(upper)
  ntl <- length(trunclower)
  ntu <- length(truncupper)
  
  # Check lengths
  if(nl!=1 & nu!=1 & nl!=nu) {
    stop("lower and upper should have equal length if both do not have length 1.")
  }
  
  if(ntl!=1 & ntu!=1 & ntl!=ntu) {
    stop("trunclower and truncupper should have equal length if both do not have length 1.")
  }
  
  
  
  # Check data types
  if(any(!is.numeric(lower) & !is.na(lower))) {
    stop("lower should consist of numerics and/or NAs.")
  }
  
  if(any(!is.numeric(upper) & !is.na(upper))) {
    stop("upper should consist of numerics and/or NAs.")
  }
  
  if(any(!is.numeric(trunclower))) {
    stop("trunclower should consist of numerics.")
  }
  
  if(any(!is.numeric(truncupper))) {
    stop("truncupper should consist of numerics.")
  }
  
  
  # Check inequalities
  if(!all(is.na(lower)) & !all(is.na(upper))) {
    if(any(lower>upper)) {
      stop("lower should be smaller than (or equal to) upper.")
    }
  }
  
  
  if(any(trunclower>truncupper)) {
    stop("trunclower should be smaller than (or equal to) truncupper.")
  }
  
  
  if(any(!is.finite(trunclower) & !is.finite(truncupper))) {
    stop("trunclower and truncupper cannot be both infinite.")
  }
  
  if(!all(is.na(lower))) {
    if(any(trunclower>lower)) {
      stop("trunclower should be smaller than (or equal to) lower.")
    }
  }
  
  
  if(!all(is.na(upper))) {
    if(any(truncupper<upper)) {
      stop("truncupper should be larger than (or equal to) cupper.")
    }
  }
  
}



#' Fit mixture of Erlangs
#' 
#' Fits univariate and multivariate mixtures of Erlang distributions to possibly
#' censored and/or truncated data. The censoring and/or truncation can be left, 
#' right or interval.
#' 
#' @param lower,upper Matrix specifying the lower and upper censoring points, observations in 
#'    rows, dimensions in coloms.
#' @param trunclower,truncupper Numeric vector specifying the lower and upper truncation 
#'    points in each dimension.
#' @param M Numeric vector of values for the tuning parameter M.
#' @param s Numeric vector of values for the tuning parameter s (the spread).
#' @param nCores Number of cores available for parallel computation.
#' @param criterium Character vector specifying information criterium to use, 
#'    either "AIC" or "BIC".
#' @param eps Numeric: covergence threshold used in the EM algorithm.
#' @param beta_tol Numeric: threshold for the mixing weights below which the 
#'    corresponding shape parameter vector is considered neglectable.
#' @param print Logical: print intermediate results, either TRUE or FALSE.
#' @param file If print is TRUE, specify file name.
#' @param max_iter Maximum number of iterations in a single EM algorithm.
#' 
#' @details 
#'    \code{MME_tune} implements the estimation procedure for univariate and 
#'    multivariate mixtures of Erlangs (MME) by repeatedly using the EM algorithm. More information 
#'    on the initialization and adjustment strategy for the shape parameter vectors based on an 
#'    information criterium (AIC and BIC implemented) can be found in Verbelen et al. (2015). 
#'    
#'    The data can be censored and/or truncated. The censoring status of a particular observation 
#'    in a certain dimension is determined as follows:
#'    \itemize{
#'        \item Uncensored: lower and upper are equal (upper is set equal to lower by default).
#'        \item Left Censored: lower is missing (NA) or equal to trunclower, but upper is present.
#'        \item Right Censored: lower is present, but upper is missing (NA) or equal to truncupper.
#'        \item Interval Censored: lower and upper are present and different.
#'    }      
#'    E.g.: lower[1, ] = c(2, NA, 4, 5) and upper[1, ] = c(2, 3, NA, 6);
#'    specifies a first four-dimensional observation having an observed event at 2 in dimension 1, 
#'    left censoring at 3 in dimension 2, 
#'    right censoring at 4 in dimension 3, 
#'    and interval censoring at [5,6] in dimension 4.
#'    
#'    The truncation status in a certain dimension is determined as follows:
#'    \itemize{
#'        \item Untruncated: trunclower equals 0 and truncupper equals Inf (default).
#'        \item Left Truncated: trunclower is a nonzero numeric and truncupper is Inf.
#'        \item Right Truncated: trunclower is 0 and truncupper is a nonzero numeric.
#'        \item Interval Truncated: trunclower and truncupper are present and different.
#'    }
#'    E.g.: trunclower = c(0, 1, 0, 1) and truncupper = c(Inf, Inf, 10, 10);
#'    specifies no truncation in dimension 1,
#'    left truncation at 1 in dimension 2, 
#'    right truncation at 10 in dimension 3, 
#'    and interval truncation at [1, 10] in dimension 4.
#'    
#'    The lower and upper truncation points are fixed, i.e. the same for all observations.
#'    
#'    For all observations and across all dimensions it must hold that 
#'    trunclower <= lower <= upper <= truncupper. 
#'   
#' @return \code{MME_tune} returns a \code{list}  with the following objects: 
#'   \describe{ 
#'   \item{best_model}{The final MME, judged to be the best according
#'   to the criterium used. Value is a list.} 
#'   \item{performances}{A matrix 
#'   summarizing the performance of the different fitted MME models (M, s, 
#'   criterium, R).} 
#'   \item{all_model}{A list containing all fitted MME models.} 
#'   }
#'   
#' @references 
#'   Verbelen, R., Gong, L., Antonio, K., Badescu, A., and Lin, X. S. 
#'   (2015). Fitting mixtures of Erlangs to censored and truncated data using
#'   the EM algorithm. ASTIN Bulletin. Accepted for publication.
#'   
#'   Verbelen, R., Antonio, K., and Claeskens, G. (2015). Multivariate mixtures
#'   of Erlangs for density estimation under censoring and truncation. Submitted
#'   for publication.
#'   
#' @examples
#'    \dontrun{
#'        data(geyser, package = "MASS")
#'        MME_tune(lower = geyser, upper = geyser, trunclower = c(0, 0),
#'        truncupper = c(Inf, Inf), M = c(5, 10, 20), s = c(seq(10, 100, 10), 200),
#'        nCores = parallel::detectCores(), criterium = "BIC", eps = 1e-03,
#'        beta_tol = 10^(-5), print=TRUE, file="log.txt", max_iter = Inf)
#'    }
#' @export
MME_tune <- function(lower, upper = lower, trunclower = rep(0, ncol(as.matrix(lower))), truncupper = rep(Inf, ncol(as.matrix(lower))), M = 10, s = 100, nCores = parallel::detectCores(), criterium = "BIC", eps = 1e-03, beta_tol = 10^(-5), print=TRUE, file="log.txt", max_iter = Inf){
  
  # Check input
  .ME_checkInput(lower=lower, upper=upper, trunclower=trunclower, truncupper=truncupper)
  
  
  tuning_parameters <- expand.grid(M, s)
  
  if(nCores==1) {
    
    if(print) writeLines(c(""), file)
    i <- 1
    
    all_model <- foreach::foreach(i = 1:nrow(tuning_parameters), 
                         .export=c("lower", "upper", "trunclower", "truncupper", "tuning_parameters", "criterium", "eps", "print", "beta_tol", "max_iter", ".MME_initial", ".MME_loglikelihood", ".MME_f", ".Estep_z", ".Estep_EX", ".Mstep_T_j", ".Mstep_T", ".theta_nlm", ".MME_em", ".MME_shape_adj", ".MME_shape_red", ".MME_fit"), 
                         .errorhandling = 'remove') %do% {
                           if(print) cat(paste("M = ", tuning_parameters[i, 1], ", s = ", tuning_parameters[i, 2], "\n"), file = file, append = TRUE)
                           suppressWarnings(.MME_fit(lower, upper, trunclower, truncupper, M = tuning_parameters[i, 1], s = tuning_parameters[i, 2], criterium, eps, FALSE))
                         }
    
    
  } else {
    
    
    #######
    # Original code
    
    cl <- parallel::makePSOCKcluster(nCores)
    parallel::clusterExport(cl, c("lower", "upper", "trunclower", "truncupper", "tuning_parameters", "criterium", "eps", "print", "beta_tol", "max_iter", ".MME_initial", ".MME_loglikelihood", ".MME_f", ".Estep_z", ".Estep_EX", ".Mstep_T_j", ".Mstep_T", ".theta_nlm", ".MME_em", ".MME_shape_adj", ".MME_shape_red", ".MME_fit"), env = environment())
    doParallel::registerDoParallel(cl)
    if(print) writeLines(c(""), file)
    i <- 1
    all_model <- foreach::foreach(i = 1:nrow(tuning_parameters), .packages='matrixStats', .errorhandling = 'remove') %dopar% {
      if(print) cat(paste("M = ", tuning_parameters[i, 1], ", s = ", tuning_parameters[i, 2], "\n"), file = file, append = TRUE)
      .MME_fit(lower, upper, trunclower, truncupper, M = tuning_parameters[i, 1], s = tuning_parameters[i, 2], criterium, eps, FALSE, beta_tol, max_iter)
    }
    parallel::stopCluster(cl) 
  }
  
  performances <- data.frame(sapply(all_model, with, M), sapply(all_model, with, s),  sapply(all_model, function(x) with(x, get(criterium))), sapply(all_model, with, R))
  colnames(performances) <- c('M', 's', criterium, 'R')
  best_index <- which.min(performances[, criterium])
  best_model <- all_model[[best_index]]
  list(best_model = best_model, performances = performances, all_model = all_model)
}
