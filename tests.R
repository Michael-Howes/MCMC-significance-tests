### The parallel method

parallel_test <-  function(x0, 
                           stat, 
                           step_forward, 
                           step_back, 
                           M = 250, 
                           L = 100) {
  
  # Calculate test statistic at observed data
  # and define vector to store MCMC test statistics
  Ts <- numeric(M + 1)
  Ts[1] <- stat(x0)
  
  # Run the Markov chain backwards to latent center
  x_current <- x0
  for (i in (1:L)) {
    x_current <- step_back(x_current)
  }
  x_star <- x_current
  
  # Repeatedly run the Markov chain forwards to sample
  # MCMC values
  for (m in (2:(M+1))) {
    x_current <- x_star
    for (i in (1:L)) {
      x_current <- step_forward(x_current)
    }
    Ts[m] <- stat(x_current)
  }
  
  # Calculate p value
  p_value <- sum(Ts >= Ts[1])/(M+1)
  
  return(list(p_value = p_value, 
              T0 = Ts[1], 
              Ts = Ts, 
              M = M, 
              L = L,
              x0 = x0,
              x_star = x_star))
}


### The serial method

serial_test <- function(x0,
                        stat, 
                        step_forward,
                        step_back,
                        M = 250, 
                        L = 100) {
  # Sample location of observed data in 
  # sequence
  m_star = sample(0:M, 1)
  
  # Calculate test statistic at observed data
  # and define vector to store MCMC test statistics
  Ts <- numeric(M+1)
  Ts[m_star + 1] <- stat(x0)
  
  # Run the Markov chain backwards from m_star and
  # calculate the test statistic on MCMC samples
  if (m_star > 0) {
    x_current <- x0
    for (m in m_star:1) {
      for (i in 1:L) {
        x_current <- step_back(x_current)
      }
      Ts[m] <- stat(x_current)
    }
  }
  
  # Run the Markov chain forwards from m_star and
  # calculate the test statistic on MCMC samples
  
  if (m_star < M) {
    x_current <- x0
    for (m in (m_star + 2):(M + 1)) {
      for (i in 1:L) {
        x_current <- step_forward(x_current)
      }
      Ts[m] <- stat(x_current)
    }
  }
  
  
  # Calculate p value
  p_value <- sum(Ts >= Ts[m_star + 1])/(M+1)
  
  return(list(p_value = p_value, 
              T_0 = Ts[m_star + 1], 
              Ts = Ts, 
              M = M, 
              L = L,
              x0 = x0,
              m_star = m_star))
}

# Standard MC Test

standard_test <- function(x0, stat, sampler, M) {
  # Calculate test statistic at observed data
  # and define vector to store MC test statistics
  Ts <- numeric(M+1)
  Ts[1] <- stat(x0)
  
  # Create Monte Carlo samples and compute test statistics
  for (i in 1:M) {
    Ts[i + 1] <- stat(sampler())
  }
  
  p_value <- sum(Ts >= Ts[1])/(M + 1)
  
  return(list(p_value = p_value, 
              T_0 = Ts[1], 
              Ts = Ts, 
              M = M, 
              x0 = x0))
}

rectangle_loop <- function(M) {
  r <- nrow(M)
  c <- ncol(M)
  
  # Sample first corner
  row_1 <- sample((1:r), 1)
  col_1 <- sample((1:c), 1)
  
  # Sample second corner
  matched_cols <- (1:c)[M[row_1,] != M[row_1, col_1]]
  if (length(matched_cols) == 0) {
    return(M)
  }
  if (length(matched_cols) == 1) {
    col_2 <- matched_cols
  } else{
    col_2 <- sample(matched_cols, 1)
  }
  
  # Sample third corner
  matched_rows <- (1:r)[M[,col_2] != M[row_1, col_2]]
  if (length(matched_rows) == 0){
    return(M)
  }
  if (length(matched_rows) == 1) {
    row_2 <- matched_rows
  } else{
    row_2 <- sample(matched_rows, 1)
  }
  
  # Check if fourth corner matches
  if (M[row_2, col_1] == M[row_1, col_1]) {
    return(M)
  }
  M[c(row_1, row_2), c(col_1, col_2)] <- 1 - M[c(row_1, row_2), c(col_1, col_2)]
  
  return(M)
}


