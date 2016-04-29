# Bivariate Normal kernel
# a is horizontal axis, b is vertical axis, r is rotation angle
gauss.elliptic <- function(x, y, mu.x = 0, mu.y = 0, s.a = 1, s.b = 1,
                           r = 0, maxrate = 1){
  rot <- zapsmall(matrix(c(cos(r), sin(r), -sin(r), cos(r)), nrow = 2))
  ab <- diag(c(s.a^2, s.b^2))
  sigma <- rot %*% ab %*% t(rot)
  siginv <- solve(sigma)
  mu <- matrix(c(mu.x, mu.y), nrow = 2)
  mat <- matrix(rbind(x, y), nrow = 2)
  return(maxrate * apply(mat, 2, function(vec){
      exp(-t(vec - mu) %*% siginv %*% (vec - mu) / 2)
    }))
}


# Parametric semivariograms
# h is the lag distance, theta = c(nugget, sill, range)
# ghat is empirical semivariogram, lags is a vector of lag distances

# Spherical
sv.sphere <- function(h, theta){
  return(theta[2] * ifelse(h<theta[3], (1.5*h/theta[3]-0.5*(h/theta[3])^3), 1))
}
# Exponential
sv.expon <- function(h, theta){
  return(theta[2] * (1 - exp(-3*h/theta[3])))
}
# Gaussian
sv.gauss <- function(h, theta){
  return(theta[2] * (1 - exp(-(3*h/theta[3])^2)))
}

# Power
# theta[2] is slope, theta[3] is exponent
sv.power <- function(h, theta){
  return(theta[2] * h^theta[3])
}

# OLS objective function
sv.ss <- function(theta, lags, ghat, model){
  theta <- ifelse(theta>0, theta, 0)
  return(sum((ghat-theta[1]-sapply(lags, model, theta = theta))^2))
}

# WLS objective function
sv.wss <- function(theta, lags, n, ghat, model){
  theta <- ifelse(theta>0, theta, 0)
  return(sum(n/2*(ghat/(theta[1]+sapply(lags, model, theta = theta))-1)^2))
}
