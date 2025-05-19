library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(Rosenbrock)

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

# Example of a generic d-dimensional function, e.g., Rosenbrock function for demonstration
my_function <- function(z, n1, n2) {
  
  xprime<-z[1]
  x<-matrix(z[-1],ncol=n1-1,byrow=TRUE)
  mu<-1
  a<-1
  b<-x*0+1
  
  x<-lapply(seq_len(ncol(x)), function(i) x[i,])
  b<-lapply(seq_len(ncol(b)), function(i) b[i,])
  
  dhybrid_own(xprime, x, a, b, mu, normalize = TRUE)
}

my_function2 <- function(z, n1, n2) {
  
  xprime<-z[1]
  x<-matrix(z[-1],ncol=n1-1,byrow=TRUE)
  mu<-1
  a<-1
  b<-x*0+1
  
  x<-lapply(seq_len(ncol(x)), function(i) x[i,])
  b<-lapply(seq_len(ncol(b)), function(i) b[i,])
  
  dhybrid(xprime, x, a, b, mu)
}

# Generate a 2D grid for each pair of variables
generate_grid <- function(range, n = 100) {
  seq(range[1], range[2], length.out = n)
}

# Evaluate the function for each pair of variables
evaluate_function_pairs <- function(f, n1 ,n2 , range, n = 100) {
  results <- list()
  
  d<- (n1-1)*n2+1
  
  # Loop over each unique pair of dimensions
  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {
      # Generate grid for variables i and j
      grid_x1 <- generate_grid(range, n)
      grid_x2 <- generate_grid(range, n)
      data <- expand.grid(x1 = grid_x1, x2 = grid_x2)
      
      # Calculate the function value for the pair (i, j) with other variables fixed at 0
      data$z <- apply(data, 1, function(row, n1, n2) {
        x <- rep(1, d)  # Fix other variables to zero
        x[i] <- row["x1"]
        x[j] <- row["x2"]
        f(x, n1, n2)
      }, n1, n2)
      
      # Label for the variable pair
      data$pair <- paste("x", i, "vs", "x", j, sep = "")
      
      # Store each result
      results[[length(results) + 1]] <- data
    }
  }
  
  # Combine all results into a single data frame
  density_df <- bind_rows(results)
  return(density_df)
}

# Parameters
n1 <- 3 
n2 <- 2
range <- c(-10, 10)  # Define the range for the grid
n <- 101  # Number of points along each dimension

# Compute function values for each pair of variables
density_df <- evaluate_function_pairs(my_function, n1, n2, range, n)
density_df2 <- evaluate_function_pairs(my_function2, n1, n2, range, n)

# Plot with ggplot
ggplot(density_df, aes(x = x1, y = x2, z = z)) +
  geom_contour() +
  facet_wrap(~ pair, scales = "free") +
  theme_minimal() +
  labs(title = "Contour Plots of Function by Variable Pair",
       x = "x1", y = "x2")

# Plot with ggplot
ggplot(density_df2, aes(x = x1, y = x2, z = z)) +
  geom_contour() +
  facet_wrap(~ pair, scales = "free") +
  theme_minimal() +
  labs(title = "Contour Plots of Function by Variable Pair",
       x = "x1", y = "x2")


a <- 1/20
mu <- 1
b <- list(c(5,5),
          c(5,5))
n <- 1e+4
z<-rhybrid(n = n,a = a,b = b,mu = mu)



# Function to create matrix of marginal density plots
create_marginal_matrix <- function(data, variables = NULL) {
  # If no variables specified, use all numeric columns
  if (is.null(variables)) {
    variables <- names(data)[sapply(data, is.numeric)]
  }
  
  # Number of variables
  n <- length(variables)
  
  # Set up the plotting matrix
  par(mfrow = c(n, n), mar = c(2, 2, 1, 1))
  
  # Create the matrix of plots
  for (i in 1:n) {
    for (j in 1:n) {
      if (i == j) {
        # Diagonal: density plots
        d <- density(data[[variables[i]]])
        plot(d, main = "", xlab = "", ylab = "")
        polygon(d, col = "lightblue", border = "blue")
      } else {
        # Off-diagonal: scatter plots
        plot(data[[variables[j]]], data[[variables[i]]], 
             pch = 20, cex = 0.6, col = "blue",
             xlab = "", ylab = "")
      }
    }
  }
  # Reset plotting parameters
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}


upperfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_density2d()+
    scale_x_continuous(limits = c(-1.5,1.5))+
    scale_y_continuous(limits = c(-1.5,1.5))
}   

lowerfun <- function(data,mapping){
  ggplot(data = data, mapping = mapping)+
    geom_point()+
    scale_x_continuous(limits = c(-1.5,1.5))+
    scale_y_continuous(limits = c(-1.5,1.5))
}  


ggpairs(as.data.frame(z),upper = list(continuous = wrap(upperfun)),
        lower = list(continuous = wrap(lowerfun)))   
