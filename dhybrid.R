#' Hybrid Rosenbrock Density Function
#' 
#' @param xprime Scalar value for the Gaussian component
#' @param x List of 2D vectors representing pairs of points
#' @param a Positive scalar, precision parameter for Gaussian component
#' @param b List of scale parameters for Rosenbrock components
#' @param mu Scalar mean parameter for Gaussian component
#' @param normalize Logical, whether to include normalizing constant
dhybrid <- function(xprime, x, a, b, mu, normalize = TRUE, log_d = FALSE) {
  # Input validation
  if (!is.numeric(xprime) || length(xprime) != 1) {
    stop("xprime must be a scalar value")
  }
  if (!is.list(x)) {
    stop("x must be a list of vectors")
  }
  if (!is.numeric(a) || length(a) != 1 || a <= 0) {
    stop("a must be a positive scalar")
  }
  if (!is.list(b) || length(b) != length(x)) {
    stop("b must be a list of same length as x")
  }
  if (!is.numeric(mu) || length(mu) != 1) {
    stop("mu must be a scalar value")
  }
  
  # Get dimensions
  n2 <- length(x)  # Number of dimension groups
  n1 <- length(x[[1]]) + 1  # Dimension within groups
  
  # Validate pairs in x
  if (!all(sapply(x, length) == n1 - 1)) {
    stop("All vectors in x must be of same length")
  }
  
  # Validate b parameters
  if (!all(sapply(b, length) == n1 -1 )) {
    stop("Each element in b must be of same length")
  }
  
  # Calculate normalizing constant
  log_const <- if (normalize) {
    ((n2 * (n1-1) + 1)/2) * log(pi) - 0.5 * log(a) - 
      0.5 * sum(unlist(lapply(X = b, FUN = log)))
  } else {
    0
  }
  
  x<-lapply(X = x, FUN = function(x) return(c(xprime,x)) )
  
  # Inner product function for Rosenbrock terms
  innerprod <- function(z, p) {
    # z should be length n1, p should be length n1-1
    
    z_sqrd<-z[1:(length(z)-1)]
    z_sqrd<-z_sqrd^2
    
    return( sum(p * (z[-1] - z_sqrd)^2) )
  }
  
  # Calculate density
  log_density <- -(a * (xprime - mu)^2 + 
                 sum(mapply(z = x, p = b, FUN = innerprod))) - log_const
  
  if(log_d == TRUE)
    return(log_density)
  else
    return(exp(log_density))
}

x<-list(c(0,0,0))
x<-c(x,list(c(0,0,0)))
mu<-1
a<-1/20
b<-list(rep(100/20,length(x[[1]])))
b<-c(b,b)
xprime<-0


dhybrid(xprime, x, a, b, mu, normalize = TRUE)
dhybrid(xprime, x, a, b, mu, normalize = TRUE, log_d=TRUE)

sqrt(a)*prod(sqrt(unlist(b)))/pi^((3*2+1)/2)*exp(-a*mu^2)
