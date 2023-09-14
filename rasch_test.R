library("difR")
library("dplyr")
library("eRm")
library("ggplot2")

source("tests.R")

set.seed(102)

# Binary matrix
data(verbal, package = "difR")
x0 <- as.matrix(select(verbal, !Gender & !Anger))
gender <- verbal$Gender

# A reversible Markov chain on the space of binary matrices
# with given row and column sums. The uniform distribution is 
# the stationary distribution of this Markov chain.

rectangle_loop <- function(A) {
  r <- nrow(A)
  c <- ncol(A)
  
  # Sample first corner
  row_1 <- sample((1:r), 1)
  col_1 <- sample((1:c), 1)
  
  # If first corner is a 1, sample second corner 
  # from available rows
  if (A[row_1, col_1] == 1){
    # Sample second corner among entries in row_1
    matched_cols <- (1:c)[A[row_1,] != A[row_1, col_1]]
    if (length(matched_cols) == 0) {
      return(A)
    }
    if (length(matched_cols) == 1) {
      col_2 <- matched_cols
    } else{
      col_2 <- sample(matched_cols, 1)
    }
    # Sample third corner among entries in col_2
    matched_rows <- (1:r)[A[,col_2] != A[row_1, col_2]]
    if (length(matched_rows) == 0){
      return(A)
    }
    if (length(matched_rows) == 1) {
      row_2 <- matched_rows
    } else{
      row_2 <- sample(matched_rows, 1)
    }
    if (A[row_2, col_1] == A[row_1, col_1]) {
      return(A)
    }
  # If first corner is a 0, sample second corner from
  # available columns
  } else {
    # Sample first corner among entries in col_1
    matched_rows <- (1:r)[A[,col_1] != A[row_1, col_1]]
    if (length(matched_rows) == 0){
      return(A)
    }
    if (length(matched_rows) == 1) {
      row_2 <- matched_rows
    }
    else {
      row_2 <- sample(matched_rows, 1)
    }
    # Sample third corner among entries in row_2 
    matched_cols <- (1:c)[A[row_2,] != A[row_2, col_1]]
    if (length(matched_cols) == 0) {
      return(A)
    }
    if (length(matched_cols) == 1) {
      col_2 <- matched_cols
    } else{
      col_2 <- sample(matched_cols, 1)
    }
    # Check if fourth corner matches
    if (A[row_1, col_2] == A[row_1, col_1]) {
      return(A)
    }
  }
  
  
  # Flipped the selected 2 by 2 sub matrix
  A[c(row_1, row_2), c(col_1, col_2)] <- 
      1 - A[c(row_1, row_2), c(col_1, col_2)]
  
  return(A)
}

# Parameters for exchangeable sampler
M <- 100
L <- 100

# Andersen's likelihood ratio test statistic
LR_stat <- function(x){
  model <- RM(x)
  return(LRtest(model, splitcr = gender)$LR)
}

result <- serial_test(x0, LR_stat, rectangle_loop, rectangle_loop, M = M, L = L)

df <- tibble(
  xs = 0:result$M,
  Ts = result$Ts
)

# Create plot
ggplot(df) + 
  geom_point(aes(xs, Ts),size = 0.5) + 
  geom_point(aes(result$m_star, result$T_0), fill = "red", size = 0.6, pch = 21) +
  geom_vline(aes(xintercept = result$m_star), color = "red") +
  labs(x = "Index", y = "Test statistic") +
  theme_bw()


print(result$p_value)
