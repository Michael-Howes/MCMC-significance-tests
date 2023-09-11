source("tests.R")
library(tidyverse)
library(latex2exp)
library(rbenchmark)
set.seed(101)


## The density pi
X <- 1:100
pi <- 0.5 * dnorm(X, 25, 6) + 0.5 * dnorm(X, 75, 6)
pi_normed <- pi/sum(pi)
CDF_pi <- cumsum(pi_normed)
density = tibble(x = X,
                 pi = pi_normed)

# 95 quantile of pi
x_star = min(X[cumsum(pi_normed) >= 0.95])
ggplot() +
  geom_line(aes(X, pi_normed),density) +
  geom_ribbon(aes(x = x_star:100, ymin = 0, ymax = pi_normed[x_star:100]),
              fill = "black") +
  theme_bw() +
  labs(y = TeX(r"($\pi(x)$)"),
       x = TeX(r"($x$)"))
ggsave("../figures/fig bimodal.pdf",
       device = "pdf",
       width = 3.3,
       height = 2.8)

## A reversible Markov chain with stationary distribution
## pi
MH_step <- function(x){
  x_prime <- x
  if (rbinom(1,1,0.5) == 1){
    x_prime <- min(x+1, 100L)
  } else {
    x_prime <- max(x-1, 1L)
  }
  if (pi[x_prime]/pi[x] < runif(1)){
    x_prime <- x
  }
  return(x_prime)
}



# An exact sampler from pi
sampler <- function() {
  X <- sample(1:100, 1, prob = pi_normed)
  return(X)
}



N_reps <- 200
x0s <- c(40,85)
N_xs = length(x0s)

L <- 100
M <- 100

p_values <- matrix(nrow = N_xs*N_reps, ncol = 6)

colnames(p_values) <- c("parallel", "serial", "standard",
                        "L", "M", "x0")

p_values[,"L"] <- L
p_values[,"M"] <- M
p_values[,"x0"] <- sapply(x0s, function(x){rep(x,N_reps)})

for (i in 1:(N_xs*N_reps)) {
  x0 <- p_values[i, "x0"]
  p_values[i, "parallel"] <- parallel_test(x0, 
                                           identity, 
                                           MH_step, 
                                           MH_step, 
                                           M, 
                                           L)$p_value
  
  
  p_values[i, "serial"] <- serial_test(x0, 
                                       identity, 
                                       MH_step, 
                                       MH_step, 
                                       M, 
                                       L)$p_value
  
  
  p_values[i, "standard"] <- standard_test(x0, 
                                           identity, 
                                           sampler,
                                           M)$p_value

}



df_ps <- as_tibble(p_values) %>% 
  pivot_longer(c("parallel", "serial", "standard"),
               names_to = "method",
               values_to = "p_value")

df_ps %>% 
  group_by(method, x0) %>% 
  summarize(power = mean(p_value <= 0.05),
            mean_p = mean(p_value),
            var_p = var(p_value)) %>% 
  arrange(x0)




