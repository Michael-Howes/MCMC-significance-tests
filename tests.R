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


# Standard Monte Carlo Test
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




