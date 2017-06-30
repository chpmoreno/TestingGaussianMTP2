library(MASS)
library(msm)
library(tidyverse)
library(reshape2)

# Naive approach
Sigma_simulation <- function(n, p) {
      Sigma <- matrix(NA, ncol = p, nrow = p)
      diag(Sigma) <- 1
      for(i in 1:p) {
            for(j in 1:p) {
                  if(is.na(Sigma[i, j]) & (i != j)) {
                        aux <- runif(1, 0, 0.5)
                        Sigma[i, j] <- aux
                        Sigma[j, i] <- aux
                  }
            }
      }
      return(Sigma)
}

Sigma_pos_def <- function(n, p) {
      while_indicator <- 1
      while(while_indicator == 1) {
            Sigma <- Sigma_simulation(n, p) 
            if(sum(eigen(Sigma)$values < 0) == 0 ){
                  while_indicator <- 0
            }
      }
      return(Sigma)
}

A_estimation <- function(p, D, I, dist = "trunc_norm", mean_trunc_norm = 0, 
                         sd_trunc_norm = 1) {
      A   <- matrix(0, p, p)
      if(dist == "trunc_norm") {
            A[p, 1:(p - 1)] <- rtnorm(p - 1, mean_trunc_norm, sd_trunc_norm, 0)
      }
      if(dist == "uniform") {
            A[p, 1:(p - 1)] <- runif(p - 1)
      }
      
      for(k in (p - 1):2) {
            for(j in 1:(k - 1)) {
                  A_cond <- 1
                  start_time <- Sys.time()
                  while(A_cond == 1) {
                        if(dist == "trunc_norm") {
                              A[k, j] <- rtnorm(1, mean_trunc_norm, sd_trunc_norm, 0)
                        }
                        if(dist == "uniform") {
                              A[k, j] <- runif(1)
                        }
                        sum_k_gt_pl <- 0
                        for(p_l in (k + 1):p) {
                              sum_k_gt_pl <- sum_k_gt_pl + ((1 / D[k, k]) * A[k, p_l] * A[k, j])
                              if(A[p_l, j] >= sum_k_gt_pl) {
                                    A_cond <- 0
                              }
                              timer <- Sys.time() - start_time
                              if(as.numeric(timer) < 60) {
                                    A[p, 1:(p - 1)] <- rtnorm(p - 1, mean(Sigma), 1, 0)
                              }
                        }
                        #if(A[k, j] >= (D[k, k] / D[k + 1, k + 1]) * A[k + 1, k] * A[k + 1, j]) {
                        #     A[p, j]
                        #     A_cond <- 0
                        #}
                  }
            }
      }
      K <- t(I - A) %*% solve(D) %*% (I - A)
      return(list(A = A, K = K))
}

n <- 1000
p <- 4
D   <- I <- diag(p) 
eps <- mvrnorm(n, rep(0,p), D)
Sigma <- Sigma_pos_def(n, p)
X   <- mvrnorm(n, rep(0, p), Sigma)


iter <- 10000
K_array <- array(NA, c(p, p, iter))
A_array <- array(NA, c(p, p, iter))
for(i in 1:iter) {
      aux_estim <- A_estimation(p, D, I, dist = "trunc_norm", mean_trunc_norm = mean(Sigma))
      K_array[ , , i] <- aux_estim$K
      A_array[ , , i] <- aux_estim$A
}

apply(A_array, 1:2, mean)
solve(Sigma)
apply(K_array, 1:2, mean)
apply(K_array, 1:2, min)
apply(K_array, 1:2, max)


K_df_names <- NULL
K_df <- NULL
for(i in 1:p) {
      for(j in 1:p) {
            if(i <=j){
                  K_df_names = c(K_df_names, paste0("K[", i, ",", j, "]"))
                  K_df <- cbind(K_df, K_array[i, j, ])
            }
      }
}
K_df <- as_data_frame(K_df)
colnames(K_df) <- K_df_names

K_df_melt <- melt(K_df)
colnames(K_df_melt)

ggplot(data = K_df_melt, aes(value, colour = variable, fill = variable, group = variable)) + 
      geom_density(alpha = 0.9) + facet_grid(variable~., scales = "free") +
      theme_dark()

