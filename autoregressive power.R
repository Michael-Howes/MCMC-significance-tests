library(tidyverse)
library(latex2exp)


# Define power function
AR_power <- function(mu, L, rho, alpha) {
  return(1 - pnorm(qnorm(1-alpha) - sqrt(1 - rho^{2*L})*mu))
}

# make plot
mu <- 2
alpha = 0.05
L <- 1:15
rho <- c(0.7, 0.85, 0.9, 0.95)

df <- expand_grid(L, rho) %>% 
  mutate(power = AR_power(mu, L, rho, alpha),
         rho_factor = factor(rho))

p <- df %>% 
  ggplot(aes(x = L, y = power, color = rho_factor)) +
  geom_line() +
  geom_segment(aes(x = min(L),
                   xend = max(L),
                   y = AR_power(mu, 1, 0, alpha), 
                   yend = AR_power(mu, 1, 0, alpha)),
               color = "black",
               linetype = "dashed", 
               linewidth = 0.5) +
  theme_bw() +
  labs(color = TeX(r"($\rho$)")) +
  theme(legend.position = "bottom")
p


ggsave("../figures/fig power.pdf",
       p,
       device = "pdf",
       width = 3.3,
       height = 3.3)
