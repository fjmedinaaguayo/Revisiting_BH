library(matrixStats)
library(mvtnorm)
library(parallel)
library(invgamma)

logTargetUnDens<-function(x,lUnDens_list,lWeights){
  
  d<-length(x)
  nDens<-length(lUnDens_list)
  lNormWeights<-lWeights-logSumExp(lWeights)
  
  lDens_eval<-sapply(1:nDens, function(i) lUnDens_list[[i]](x,d=d))
  
  val<-logSumExp(lDens_eval+lNormWeights)
  
  return(val)
}

ld0<-function(x,type="lUnDen",d){
  
  mu<-rep(0,d)
  diagVec<-rep(5,d)
  
  Sigma<-diag(diagVec)
  Sigma_inv<-diag(1/diagVec)
  
  if(type=="lZ")
    return(d/2*log(2*pi)+0.5*log(abs(det(Sigma))))
  if(type=="draw")
    return(rmvnorm(1,mu,Sigma))
  else if(type=="lUnDen")
    return(-0.5*sum((x-mu)%*%Sigma_inv%*%(x-mu)))
}

ld1<-function(x,type="lUnDen",d){
  
  Transf<-function(x,b){
    
    x[2]<-x[2]+b*x[1]^2-100*b
    
    return(x)
  }
  
  mu<-rep(0,d)
  diagVec<-rep(1,d)
  diagVec[1]<-diagVec[1]*100
  Sigma<-diag(diagVec)
  Sigma_inv<-diag(1/diagVec)
  
  if(type=="lZ")
    return(d/2*log(2*pi)+0.5*log(abs(det(Sigma))))
  
  b<-0.03
  x<-Transf(x,b)
  
  if(type=="draw")
    return(NA)
  else if(type=="lUnDen")
    return(-0.5*sum((x-mu)%*%Sigma_inv%*%(x-mu)))
}

ld2<-function(x,type="lUnDen",d){
  
  mu<-rep(-10,d)
  Sigma<-diag(rep(1,d))
  Sigma[1,2]<-2
  Sigma[2,1]<-2
  Sigma[2,2]<-5
  Sigma_inv<-solve(Sigma)
  
  if(type=="draw")
    return(rmvnorm(1,mu,Sigma))
  else if(type=="lZ")
    return(d/2*log(2*pi)+0.5*log(abs(det(Sigma))))
  else if(type=="lUnDen")
    return(-0.5*sum((x-mu)%*%Sigma_inv%*%(x-mu)))
}

ld3<-function(x,type="lUnDen",d){
  
  nu<-3
  mu<-rep(5,d)
  Sigma<-diag(rep(1,d))
  Sigma_inv<-solve(Sigma)
  
  if(type=="draw")
    return(rmvt(1,sigma=Sigma,df=nu,delta=mu))
  else if(type=="lZ")
    return(lgamma(nu/2)-lgamma(0.5*(nu+d))+d/2*log(nu*pi)+0.5*log(abs(det(Sigma))))
  else if(type=="lUnDen")
    return(-0.5*(nu+d)*log(1+1/nu*sum((x-mu)%*%Sigma_inv%*%(x-mu))))
}


ld4<-function(x,type="lUnDen",d){
  
  nu<-3
  mu<-rep(5,d)
  Sigma<-diag(rep(1,d))
  Sigma_inv<-solve(Sigma)
  
  if(type=="draw")
    return(rmvt(1,sigma=Sigma,df=nu,delta=mu))
  else if(type=="lZ")
    return(lgamma(nu/2)-lgamma(0.5*(nu+d))+d/2*log(nu*pi)+0.5*log(abs(det(Sigma))))
  else if(type=="lUnDen")
    return(-0.5*(nu+d)*log(1+1/nu*sum((x-mu)%*%Sigma_inv%*%(x-mu))))
}

logSumExp_by_rows <- function(x) {
  # Handle case where entire row is -Inf
  all_neg_inf <- apply(x, 1, function(row) all(is.infinite(row) & row < 0))
  
  # For each row, find the maximum value in that row
  row_max <- apply(x, 1, max)
  
  # If a row's max is -Inf, that means all elements are -Inf
  # logSumExp of all -Inf should be -Inf
  if(any(is.infinite(row_max) & row_max < 0)) {
    browser()
    result <- row_max  # Will be -Inf for rows where all elements are -Inf
  } else {
    # Proceed with normal computation for other rows
    # Subtract the row maximum from each element (for numerical stability)
    row_max_matrix <- matrix(rep(row_max, ncol(x)), nrow = nrow(x), byrow = FALSE)
    stable_exp_sum <- rowSums(exp(x - row_max_matrix))
    
    # Compute logsumexp by adding the log of the summed exponentials
    result <- row_max + log(stable_exp_sum)
  }
  
  # Set result to -Inf for rows where all values were -Inf
  result[all_neg_inf] <- -Inf
  
  if(any(is.na(result)))
    return(result)
  
  if(any(result==Inf) || any(is.nan(result)))
    browser()
  
  return(result)
}


ld_mixture_model<-function(x,type='lUnDen',d){
  
  k<-d/2
  
  mu_i<-x[1:k]
  log_sigma_i<-x[(k+1):(2*k)]
  sigma_i_sq<-exp(2*log_sigma_i)
  
  if(any(diff(mu_i)<0)){
    #browser()
    return(-Inf)
  }
  
  log_prior_mu <- sum(dnorm(mu_i, mean = mu_0, sd = tau, log = TRUE))
  log_prior_sigma_sq <- sum(dinvgamma(sigma_i_sq, shape = a, rate = b, log = TRUE))
  
  log_likelihoods <- sapply(1:k, function(i) {
    dnorm(y, mean = mu_i[i], sd = sqrt(sigma_i_sq[i]), log = TRUE)-log(k)
  })
  
  log_likelihood <- sum(logSumExp_by_rows(log_likelihoods))
  
  res<-log_prior_mu + log_prior_sigma_sq + log_likelihood
  
  if(is.null(res))
    browser()
  
  if(is.nan(res))
    browser()
  
  return(res)
  
}

ld_rosenbrock<-function(z,type='lUnDen',d){
  
  dhybrid_own <- function(xprime, x, a, b, mu, normalize = TRUE, log_d = FALSE) {
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
  
  n1<-3
  n2 <- (d-1)/(n1-1) + 1
  mu<-1
  #a<-1/20
  a<-1
  
  if(type=="lUnDen"){
    
    
    xprime<-z[1]
    x<-matrix(z[-1],ncol=n1-1,byrow=TRUE)
    #b<-x*0+100/20
    b<-x*0+1
    
    x<-lapply(seq_len(nrow(x)), function(i) x[i,])
    b<-lapply(seq_len(nrow(b)), function(i) b[i,])
    res <- dhybrid_own(xprime, x, a, b, mu, normalize = FALSE, log_d = TRUE)
    
    return(res) 
  }
  else if(type=="lZ"){
    
    xprime<-z[1]*0
    x<-matrix(z[-1]*0,ncol=n1-1,byrow=TRUE)
    #b<-x*0+100/20
    b<-x*0+1
    
    x<-lapply(seq_len(nrow(x)), function(i) x[i,])
    b<-lapply(seq_len(nrow(b)), function(i) b[i,])
    res <- dhybrid_own(xprime, x, a, b, mu, normalize = FALSE, log_d = TRUE) - 
      dhybrid_own(xprime, x, a, b, mu, normalize = TRUE, log_d = TRUE)
    
    return(res)
  }

}


ld_my_rosenbrock<-function(z,type='lUnDen',d){
  
  a<-0
  b<-100
  c<-0.1
  
  my_rosenbrock<-function(x,a=0,b=100,c=0.1,log=FALSE){
    
    d<-length(x)
    y<-x^2
    
    term1<-b*sum((x[2:d]-y[1:(d-1)])^2)
    term2<-sum((a-x[1:(d-1)])^2)
    
    res<- -0.5/c*(term1+term2)
    
    if(log==TRUE)
      return(res)
    else
      return(exp(res))
  }
  
  if(type=="lUnDen"){
    
    res <- my_rosenbrock(z, log = TRUE)
    
    return(res) 
  }
  else if(type=="lZ"){
    
    res <- 0.5*d*log(2*pi*c) -0.5*log(b) - 0.5*(d-2)*log(b+1)
    
    return(res)
  }
  
}


