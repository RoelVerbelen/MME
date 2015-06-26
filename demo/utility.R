devAskNewPage(ask = FALSE)
## Univariate mixture of Erlangs density function

ME_density <- function(x, theta, shape, alpha, trunclower = 0, truncupper = Inf, log = FALSE){
  f <- outer(x, shape, dgamma, scale = theta)
  d <- rowSums(t(t(f)*alpha))
  if(!(trunclower==0 & truncupper==Inf)){
    d <- d / (ME_cdf(truncupper, theta, shape, alpha) - ME_cdf(trunclower, theta, shape, alpha)) * ((trunclower <= x) & (x <= truncupper))
  }
  if(log){
    d <- log(d)
  }
  d
}

## Multivariate mixture of Erlangs density function
# x: Vector or matrix of observations. If x is a matrix, each row is taken to be an observation.

MME_density <- function(x, theta, shape, alpha, trunclower=rep(0, ncol(shape)), truncupper=rep(Inf, ncol(shape))){
  if(is.vector(x)){
    x <- matrix(x, ncol = ncol(shape))
  }
  f <- matrix(0, nrow = nrow(x), ncol = nrow(shape))
  for(j in 1:nrow(shape)){
    f[,j] = rowProds(t(dgamma(t(x), shape=shape[j,], scale=theta)))
  }
  d <- rowSums(t(t(f)*alpha))
  if(any(c(trunclower!=0, truncupper!=Inf))){
     t_probabilities <- colProds(matrix( pgamma(truncupper, shape=t(shape), scale=theta) - pgamma(trunclower, shape=t(shape), scale=theta) , dim(shape)[2], dim(shape)[1]))
     d <- d / sum(alpha *t_probabilities) * apply(x, 1, function(x, trunclower, truncupper){all(c(trunclower <= x, x <= truncupper))}, trunclower, truncupper)
  }
  d
}

## Univariate mixture of Erlangs cumulative distribution function

ME_cdf <- function(x, theta, shape, alpha, trunclower = 0, truncupper = Inf, lower.tail = TRUE, log.p = FALSE){
  cdf <- outer(x, shape, pgamma, scale=theta)
  p <- rowSums(t(t(cdf)*alpha))
  if(!(trunclower==0 & truncupper==Inf)){
    l <- ME_cdf(trunclower, theta, shape, alpha)
    u <- ME_cdf(truncupper, theta, shape, alpha)
    p <- ((p - l) / (u - l)) ^ {(x <= truncupper)} * (trunclower <= x)
  }
  if(!lower.tail){
    p <- 1 - p
  }
  if(log.p){
    p <- log(p)
  }
  p
}

## Multivariate mixture of Erlangs cumulative distribution function
# lower/upper: Vector or matrix of observations. If lower/upper is a matrix, each row is taken to be an observation.

MME_cdf <- function(lower, upper, theta, shape, alpha, trunclower=rep(0, ncol(shape)), truncupper=rep(Inf, ncol(shape))){
  if( !all( t(t(lower) >= trunclower) & t(t(lower) <= truncupper) & t(t(upper) >= trunclower) & t(t(upper) <= truncupper) ) ){
    stop("upper and lower must lie between truncation bounds")
  }
  if(is.vector(lower) | is.vector(upper)){
    lower <- matrix(lower, ncol = ncol(shape))
    upper <- matrix(upper, ncol = ncol(shape))
  }
  f <- matrix(0, nrow = nrow(lower), ncol = nrow(shape))
  for(j in 1:nrow(shape)){
    f[,j] = rowProds(ifelse(lower == upper, t(dgamma(t(lower), shape=shape[j,], scale=theta)), t(pgamma(t(upper), shape=shape[j,], scale=theta)) - t(pgamma(t(lower), shape=shape[j,], scale=theta))))
  }
  p <- rowSums(t(t(f)*alpha))
  if(any(c(trunclower!=0, truncupper!=Inf))){
     t_probabilities <- colProds(matrix( pgamma(truncupper, shape=t(shape), scale=theta) - pgamma(trunclower, shape=t(shape), scale=theta) , dim(shape)[2], dim(shape)[1]))
     p <- p / sum(alpha *t_probabilities)
  }
  p
}

## Noncentral moments of order k (k can be vector)

ME_moment <- function(k, theta, shape, alpha){
  alpha_factorials <- t(exp(t(lgamma(outer(k,shape,'+'))) - lgamma(shape)) * alpha)
  rowSums(alpha_factorials)*theta^k
}

## Value-at-Risk (VaR) or quantile function

ME_VaR <- function(p, theta, shape, alpha, trunclower = 0, truncupper = Inf, interval = if(trunclower == 0 & truncupper == Inf){c(qgamma(p, shape = min(shape), scale = theta), qgamma(p, shape = max(shape), scale = theta))}else{c(trunclower, min(truncupper, trunclower + qgamma(p, shape = max(shape), scale = theta)))}, start = qgamma(p, shape = shape[which.max(alpha)], scale = theta)){
  if(p==1){
   return(Inf)
  }
  if(length(shape) == 1 & trunclower == 0 & truncupper == Inf){
    VaR <- qgamma(p, shape = shape, scale = theta)
  }
  else{
    objective <- function(x){return(10000000*(ME_cdf(x, theta, shape, alpha, trunclower, truncupper)-p)^2)}
    VaR_nlm <- nlm(f = objective, p = start)
    VaR_optimize <- optimize(f = objective, interval = interval)
    VaR <- ifelse(VaR_nlm$minimum < VaR_optimize$objective, VaR_nlm$estimate, VaR_optimize$minimum)
  }
  if(objective(VaR)>1e-06){ # in case optimization fails, retry with more different starting values
    alpha <- alpha[order(shape)]
    shape <- shape[order(shape)]
    VaR_nlm <-  vector("list", length(shape))
    VaR_optimize <-  vector("list", length(shape))
    interval <- c(0, qgamma(p, shape, scale = theta))
    for(i in 1:length(shape)){
      VaR_nlm[[i]] <- nlm(f = objective, p = qgamma(p, shape = shape[i], scale = theta))
      VaR_optimize[[i]] <- optimize(f = objective, interval = interval[c(i, i+1)])
    }
    VaR_nlm <- sapply(VaR_nlm, with, estimate)[which.min(sapply(VaR_nlm, with, minimum))]
    VaR_optimize <- sapply(VaR_optimize, with, minimum)[which.min(sapply(VaR_optimize, with, objective))]
    VaR <- ifelse(objective(VaR_nlm) < objective(VaR_optimize), VaR_nlm, VaR_optimize)
  }
  VaR
}

ME_VaR <- Vectorize(ME_VaR, vectorize.args = c("p", "start"))

## Tail-Value-at-Risk (TVaR)

# using a Newton-type algorithm

ME_TVaR <- function(p, theta, shape, alpha, interval = c(qgamma(p, shape = min(shape), scale = theta), qgamma(p, shape = max(shape), scale = theta)), start = qgamma(p, shape = shape[which.max(alpha)], scale = theta)){
  VaR <- ME_VaR(p, theta, shape, alpha, interval = interval, start = start)
  M <- max(shape)
  shapes <- seq(1, M)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  A <- rev(cumsum(rev(alphas)))
  # note: A_n = A[n+1]
  AA <- rev(cumsum(rev(A)))
  # note: AA_n = AA[n+1]
  VaR + theta^2/(1-p)*sum(AA*dgamma(VaR, shapes, scale=theta))
}

ME_TVaR <- Vectorize(ME_TVaR, vectorize.args = c("p", "start"))

## Expected Shortfall (ESF)

ME_ESF <- function(R, theta, shape, alpha){
  M <- max(shape)
  shapes <- seq(1, M)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  A <- rev(cumsum(rev(alphas)))
  # note: A_n = A[n+1]
  AA <- rev(cumsum(rev(A)))
  # note: AA_n = AA[n+1]
  theta^2*sum(AA*dgamma(R, shapes, scale=theta))
}

ME_ESF <- Vectorize(ME_ESF, vectorize.args = c("R"))

## Excess loss or residual lifetime distribution

ME_excess_loss <- function (R, theta, shape, alpha) {
  # add zero-weight components to mixture to ease further compuation
  M <- max(shape)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  # note: A_n = A[n+1]
  A <- rev(cumsum(rev(alphas)))
  shape_el <- 1:M
  alpha_el <- rep(0, M)
  f <- dgamma(R, shape_el, scale = theta)
  denom <- sum(A*f)
  for(i in 1:M){
    alpha_el[i] <- sum(alphas[i:M]*f[1:(M-i+1)]) / denom
  }
  list(theta = theta, shape = shape_el, alpha = alpha_el)
}

## Random generation of univariate mixture of Erlangs

ME_random <- function(n, theta, shape, alpha){
  rgamma(n, shape=sample(shape, size=n, replace=TRUE, prob=alpha), scale=theta)
}

## Random generation of multivariate mixture of Erlangs

MME_random <- function(n, theta, shape, alpha){
  R <- dim(shape)[1]
  d <- dim(shape)[2]
  shape_sample = sample(1:R, size=n, replace=TRUE, prob=alpha)
  matrix(rgamma(n*d, shape=shape[shape_sample,], scale=theta), ncol=d)
}

## Obtain the parameters of the marginal univariate or multivariate Erlang mixture
# j is the vector indicating the marginal dimension numbers

MME_marg <- function(j, theta, shape, alpha){
  shape_marg <- unique(shape[, j, drop=FALSE])
  alpha_marg <- numeric(nrow(shape_marg))
  for(i in 1:nrow(shape_marg)){
    indices_logical <- rowAlls(shape[, j, drop=FALSE] == matrix(shape_marg[i, ], nrow = nrow(shape), ncol=length(j), byrow=TRUE))
    alpha_marg[i] <- sum(alpha[indices_logical])
  }
  list(theta = theta, shape = shape_marg[, , drop = ifelse(length(j)==1, TRUE, FALSE)], alpha = alpha_marg)
}

## Obtain the parameters of the conditonal marginal univariate or multivariate Erlang mixture
# j is the vector indicating the marginal dimension numbers
# cond_values is the vector containing the conditioning values in the other dimensions

MME_cond <- function(j, theta, shape, alpha, cond_values){
  dens_cond <- numeric(length(alpha))
  for(i in 1:length(alpha)){
    dens_cond[i] <- MME_density(cond_values, theta, shape[i, -j, drop = FALSE], alpha)
  }
  alpha_c <- alpha * dens_cond / sum(alpha * dens_cond)
  shape_cond <- unique(shape[, j, drop=FALSE])
  alpha_cond<- numeric(nrow(shape_cond))
  for(i in 1:nrow(shape_cond)){
    indices_logical <- rowAlls(shape[, j, drop=FALSE] == matrix(shape_cond[i, ], nrow = nrow(shape), ncol=length(j), byrow=TRUE))
    alpha_cond[i] <- sum(alpha_c[indices_logical])
  }
  list(theta = theta, shape = shape_cond[, , drop=TRUE], alpha = alpha_cond)
}

# Obtain the parameters of the sum of the marginal univariate Erlang mixtures

MME_sum <- function(theta, shape, alpha){
  sum_shape = sort(unique(rowSums(shape)))
  sum_alpha = sapply(sum_shape, function(x, ...){sum(alpha[rowSums(shape)==x])}, alpha = alpha)
  list(theta = theta, shape = sum_shape, alpha = sum_alpha)
}

## Plot density of mixture of erlangs

ME_plot_density <- function(theta, shape, alpha, trunclower = 0, truncupper = Inf, xlim = c(0, 10), ...) {
  x <-  seq(xlim[1], xlim[2], by=(xlim[2]-xlim[1])/10000)
  y <- ME_density(x, theta, shape, alpha, trunclower, truncupper)
  plot(x, y, type="l", xlim = xlim, ...)
}

## Plot cdf of mixture of erlangs

ME_plot_cdf <- function(theta, shape, alpha, trunclower = 0, truncupper = Inf, lower.tail = TRUE, xlim = c(0, 10), ...) {
  x <-  seq(xlim[1], xlim[2], by=(xlim[2]-xlim[1])/10000)
  y <- ME_cdf(x, theta, shape, alpha, trunclower, truncupper, lower.tail)
  plot(x, y, type="l", xlim = xlim, ...)
}

## Plot density of mixture of erlangs and histogram of data

ME_plot_data <- function(data, shape, alpha, theta, trunclower = 0, truncupper = Inf, nbins = 50, xlim = c(max(c(min(data)-0.1*(max(data)-min(data)),0)), max(data)+0.1*(max(data)-min(data))), xlab = "", legend = TRUE, lwd = 2, ...) {
  x <-  seq(xlim[1], xlim[2], by=(xlim[2]-xlim[1])/10000)
  y <- ME_density(x, theta, shape, alpha, trunclower, truncupper)
  truehist(data, h = (xlim[2]-xlim[1])/nbins, x0 = trunclower, xlim = xlim, xlab = xlab, ... )
  lines(x, y, col = "red", lwd = lwd)
  if(legend){
    legend('topright', legend = if(trunclower==0 & truncupper==Inf) c("Fitted Density Function", "Observed Relative Frequency") else c("Fitted Truncated Density Function", "Observed Relative Frequency"), col = c("red","cyan"), pch=c(NA,15), pt.cex=2, lty = c(19,NA), lwd=c(lwd,NA))
  }
}

## Plot density of mixture of erlangs, histogram of simulated data and true density

ME_plot_sim_data <- function(data, dens, dens_param, shape, alpha, theta, trunclower = 0, truncupper = Inf, nbins = 50, xlim = c(max(c(min(data)-0.1*(max(data)-min(data)),0)), max(data)+0.1*(max(data)-min(data))), xlab = "", legend = TRUE, lwd = 2, ...) {
  x <-  seq(xlim[1], xlim[2], by=(xlim[2]-xlim[1])/10000)
  y <- ME_density(x, theta, shape, alpha, trunclower, truncupper)
  yy <- do.call(dens, c(list(x = x), dens_param))
  truehist(data, h = (xlim[2]-xlim[1])/nbins, x0 = trunclower, xlim = xlim, xlab = xlab, ... )
  lines(x, y, col = "red", lwd = lwd)
  lines(x, yy, col = "blue", lwd = lwd, lty=2)
  if(legend){
    legend('topright', legend = if(trunclower==0 & truncupper==Inf) c("Fitted Density Function", "True Density Function", "Observed Relative Frequency") else c("Fitted Truncated Density Function", "True Density Function", "Observed Relative Frequency"), col = c("red", "blue", "cyan"), pch=c(NA, NA,15), pt.cex=2, lty = c(1,2,NA), lwd=c(lwd, lwd,NA))
  }
}

## Stop loss moments of order delta (delta doesn't have to be an integer), given lower and upper truncation bounds
## Equal to ESF for delta = 1, trunclower = 0, truncupper = Inf

ME_SLM <- function(R, delta = 1, theta, shape, alpha, trunclower = 0, truncupper = Inf){
  M <- max(shape)
  shapes <- seq(1, M)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  coeff <- rep(0, M)
  for(n in 1:M){
    coeff[n] <- sum( alphas[0:(M-n)+n] * exp(lgamma(delta+0:(M-n)+1) - lgamma(0:(M-n)+1)) * (pgamma(truncupper-R, shape = delta+0:(M-n)+1, scale = theta) - pgamma(max(trunclower, R)-R, shape = delta+0:(M-n)+1, scale = theta)) )
  }
  theta^(delta+1) / (ME_cdf(truncupper, theta, shape, alpha) - ME_cdf(trunclower, theta, shape, alpha)) * sum( coeff * dgamma(R, shapes, scale=theta) )
}

## Excess-of-loss reinsurance premium: C xs R (Retention R, Cover C, Limit L = R+C)

ME_XL <- function(R, C, theta, shape, alpha){
  M <- max(shape)
  shapes <- seq(1, M)
  alphas <- rep(0, M)
  alphas[shape] <- alpha
  coeff <- rep(0, M)
  for(n in 1:M){
    coeff[n] <- sum( alphas[0:(M-n)+n] * (0:(M-n)+1) * pgamma(C, shape = 0:(M-n)+2, scale = theta) )
  }
  if(C == Inf){
    XL <- theta^2 * sum( coeff * dgamma(R, shapes, scale=theta) )
  }else{
    XL <- theta^2 * sum( coeff * dgamma(R, shapes, scale=theta) ) +  C * mix.erlang.cdf(R+C, theta, shape, alpha, lower.tail = FALSE)
  }
  XL
}

## Joint moment. n = (n_1, ..., n_d) vector

MME_moment <- function(n, theta, shape, alpha){
    sum(alpha*rowProds(exp(lgamma(t(t(shape)+n)) - lgamma(shape))))*theta^(sum(n))
}

## Kendall's tau and Spearman's rho for bivariate mixtures of Erlangs

# Auxiliary function
MME_Q <- function(i, j, shape, alpha){
  sum(alpha[colAlls(t(shape) > c(i, j))])
}

MME_cor <- function(theta, shape, alpha, method =  "kendall"){
  maxs <- colMaxs(shape)
  Q <- matrix(0, maxs[1], maxs[2])
  # formula (3.4) in Lee and Lin (2012)
  # note Q_ij = Q[i+1,j+1]
  for(i in 0:(maxs[1]-1)){
    for(j in 0:(maxs[2]-1)){
      Q[i+1, j+1] <- MME_Q(i, j, shape, alpha)
    }
  }
  if(method == "kendall"){
    summand <- rep(0, length(alpha))
    for(r in 1:length(alpha)){
      k <- shape[r, 1]
      l <- shape[r, 2]
      for(i in 0:(maxs[1]-1)){
        for(j in 0:(maxs[2]-1)){
          summand[r] <- summand[r] + exp(lchoose(i+k-1, i) + lchoose(j+l-1, j) + log(Q[i+1, j+1]) -(i+j+k+l)*log(2))
        }
      }
    }
    intgrl <- sum(alpha*summand)
    return(4 * intgrl - 1)
  }
  if(method == "spearman"){
    marg1 <- MME_marg(1, theta, shape, alpha)
    marg2 <- MME_marg(2, theta, shape, alpha)
    indices <- expand.grid(1:length(marg1$shape), 1:length(marg2$shape))
    R <- nrow(indices)
    summand <- 0
    intgrl <- 0
    for(r in 1:R){
      k <- marg1$shape[indices[r, 1]]
      l <- marg2$shape[indices[r, 2]]
      for(i in 0:(maxs[1]-1)){
        for(j in 0:(maxs[2]-1)){
          summand <- summand + exp(lchoose(i+k-1, i) + lchoose(j+l-1, j) + log(Q[i+1, j+1]) -(i+j+k+l)*log(2))
        }
      }
      intgrl <- intgrl + summand * marg1$alpha[indices[r, 1]] * marg2$alpha[indices[r, 2]]
      summand <- 0
    }
    return(12 * intgrl - 3)
  }
  else{
    stop("Specify method as either 'kendall' or 'spearman'")
  }
}

## Multivariate excess loss

MME_excess_loss <- function(D, theta, shape, alpha){
  d <- ncol(shape)
  # add zero-weight components to mixture to ease further compuation
  M <- colMaxs(shape)
  r <- vector("list", d)
  for(j in 1:d){
    r[[j]] <- 1:M[j]
  }
  shapes <- as.matrix(expand.grid(r))
  alphas <- rep(0, nrow(shapes))
  for(i in 1:length(alpha)){
    alphas[colAlls(t(shapes) == shape[i, ])] <- alpha[i]
  }
  alpha_el <- rep(0, nrow(shapes))
  for(i in 1:length(alpha_el)){
    suppressWarnings(alpha_el[i] <- sum( ( alphas * theta^d * colProds( dgamma(D, t(shapes) - shapes[i, ] + 1, scale = theta) ) ) [ colAlls(t(shapes) >= shapes[i, ]) ] ) )
  }
  alpha_el <- alpha_el / MME_cdf(D, rep(Inf, length(D)), theta, shape, alpha)
  list(theta = theta, shape = shapes, alpha = alpha_el)
}
