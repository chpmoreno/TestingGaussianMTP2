library(mvtnorm)
library(tidyverse)
library(reshape2)
library(matrixcalc)
library(gRbase)
library(mixAK)
library(xtable)

source("~/Dropbox/Documents/BGSE/Thesis/GM/Simulations/Test_Code/Gibbs_Truncation.R")
source("~/Dropbox/Documents/BGSE/Thesis/GM/Simulations/Test_Code/MH.R")
source("~/Dropbox/Documents/BGSE/Thesis/GM/Simulations/Test_Code/Utils.R")

# Initial common parameters
ngibbs <- 10000
nMH    <- 10000
# Initial toy Examples ####
# p = 3, n = 1000
p = 10
n = 1000

# M-Matrix (K_ij = -1/p)
K_Mm <- matrix(0, ncol = p, nrow = p)
diag(K_Mm) <- 1
for(i in 1:p) {
      for(j in 1:p) {
            if(i != j) K_Mm[i, j] <- -1/p
      }
}
Sigma_Mm <- solve(K_Mm)

X_Mm <- rmvnorm(n, rep(0,p), Sigma_Mm)

gibbs_trunc_X_Mm <- gibbs_trunc(K_Mm, X_Mm, ngibbs, burn_rate = 0.1)
gibbs_trunc_X_Mm$dens_plot_1
gibbs_trunc_X_Mm$dens_plot_2
gibbs_trunc_X_Mm$Gibbs_summary$mean
gibbs_trunc_X_Mm$is_m_matrix_prop
gibbs_trunc_X_Mm$n_hold

MH_X_Mm <- MH(K_Mm, X_Mm, nMH, burn_rate = 0.1)
MH_X_Mm$dens_plot_1
MH_X_Mm$dens_plot_2
MH_X_Mm$Gibbs_summary$mean
MH_X_Mm$is_m_matrix_prop

latex_gibbs_Mm <- xtable(gibbs_trunc_X_Mm$Gibbs_summary$mean, 
                         align=rep("", ncol(gibbs_trunc_X_Mm$Gibbs_summary$mean) + 1))
print(latex_gibbs_Mm, floating=FALSE,
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

latex_MH_Mm <- xtable(MH_X_Mm$Gibbs_summary$mean, 
                      align=rep("", ncol(MH_X_Mm$Gibbs_summary$mean) + 1))
print(latex_MH_Mm, floating=FALSE,
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

# Non-M-Matrix (K_ij = 1/p)
K_nMm <- matrix(0, ncol = p, nrow = p)
diag(K_nMm) <- 1
for(i in 1:p) {
      for(j in 1:p) {
            if(i != j) K_nMm[i, j] <- 1/p
      }
}
Sigma_nMm <- solve(K_nMm)

X_nMm <- rmvnorm(n, rep(0,p), Sigma_nMm)

gibbs_trunc_X_nMm <- gibbs_trunc(K_nMm, X_nMm, ngibbs, burn_rate = 0.1)
gibbs_trunc_X_nMm$dens_plot_1
gibbs_trunc_X_nMm$dens_plot_2
gibbs_trunc_X_nMm$Gibbs_summary$mean
gibbs_trunc_X_nMm$is_m_matrix_prop
gibbs_trunc_X_nMm$n_hold

MH_X_nMm <- MH(K_nMm, X_nMm, nMH, burn_rate = 0.1)
MH_X_nMm$dens_plot_1
MH_X_nMm$dens_plot_2
MH_X_nMm$Gibbs_summary$mean
MH_X_nMm$is_m_matrix_prop

latex_gibbs_nMm <- xtable(gibbs_trunc_X_nMm$Gibbs_summary$mean, 
                         align=rep("", ncol(gibbs_trunc_X_nMm$Gibbs_summary$mean) + 1))
print(latex_gibbs_nMm, floating=FALSE,
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

latex_MH_nMm <- xtable(MH_X_nMm$Gibbs_summary$mean, 
                      align=rep("", ncol(MH_X_nMm$Gibbs_summary$mean) + 1))
print(latex_MH_nMm, floating=FALSE,
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

# M-Matrix (K_ij = -1/p, K_1p = K_p1 = 0) 
K_Mm <- matrix(0, ncol = p, nrow = p)
diag(K_Mm) <- 1
for(i in 1:p) {
      for(j in 1:p) {
            if(i != j) K_Mm[i, j] <- -1/p
      }
}
K_Mm[1, p] <- 0
K_Mm[p, 1] <- 0

Sigma_Mm <- solve(K_Mm)

X_Mm <- rmvnorm(n, rep(0,p), Sigma_Mm)

gibbs_trunc_X_Mm <- gibbs_trunc(K_Mm, X_Mm, ngibbs, burn_rate = 0.1)
gibbs_trunc_X_Mm$dens_plot_1
gibbs_trunc_X_Mm$dens_plot_2
gibbs_trunc_X_Mm$Gibbs_summary$mean
gibbs_trunc_X_Mm$is_m_matrix_prop
gibbs_trunc_X_Mm$n_hold

MH_X_Mm <- MH(K_Mm, X_Mm, nMH, burn_rate = 0.1)
MH_X_Mm$dens_plot_1
MH_X_Mm$dens_plot_2
MH_X_Mm$Gibbs_summary$mean
MH_X_Mm$is_m_matrix_prop


# Non-M-Matrix (K_ij = 1/p, K_1p = K_p1 = 0)
K_nMm <- matrix(0, ncol = p, nrow = p)
diag(K_nMm) <- 1
for(i in 1:p) {
      for(j in 1:p) {
            if(i != j) K_nMm[i, j] <- 1/p
      }
}
K_nMm[1, p] <- 0
K_nMm[p, 1] <- 0

Sigma_nMm <- solve(K_nMm)

X_nMm <- rmvnorm(n, rep(0,p), Sigma_nMm)

gibbs_trunc_X_nMm <- gibbs_trunc(K_nMm, X_nMm, ngibbs, burn_rate = 0.1)
gibbs_trunc_X_nMm$dens_plot_1
gibbs_trunc_X_nMm$dens_plot_2
gibbs_trunc_X_nMm$Gibbs_summary$mean
gibbs_trunc_X_nMm$is_m_matrix_prop
gibbs_trunc_X_nMm$n_hold

MH_X_nMm <- MH(K_nMm, X_nMm, nMH, burn_rate = 0.1)
MH_X_nMm$dens_plot_1
MH_X_nMm$dens_plot_2
MH_X_nMm$Gibbs_summary$mean
MH_X_nMm$is_m_matrix_prop

# Time running:
p_vector <- c(3, 4, 5, 10, 20)
time_p_vector_gibbs <- NULL
for(t in p_vector) {
      print(t)
      time_p_ini <- Sys.time()
      p = t
      n = 1000

      # M-Matrix (K_ij = -1/p)
      K_Mm <- matrix(0, ncol = p, nrow = p)
      diag(K_Mm) <- 1
      for(i in 1:p) {
            for(j in 1:p) {
                  if(i != j) K_Mm[i, j] <- -1/p
            }
      }
      Sigma_Mm <- solve(K_Mm)

      X_Mm <- rmvnorm(n, rep(0,p), Sigma_Mm)
      gibbs_trunc_X_Mm <- gibbs_trunc(K_Mm, X_Mm, ngibbs, burn_rate = 0.1, plot_dens = FALSE, statistics = FALSE)
      time_p_fin <- Sys.time()
      time_p_vector_gibbs <- c(time_p_vector_gibbs, difftime(time_p_fin, time_p_ini, units = "mins"))
}

p_vector <- c(3, 4, 5, 10, 20)
time_p_vector_MH <- NULL
for(t in p_vector) {
      print(t)
      time_p_ini <- Sys.time()
      p = t
      n = 1000

      # M-Matrix (K_ij = -1/p)
      K_Mm <- matrix(0, ncol = p, nrow = p)
      diag(K_Mm) <- 1
      for(i in 1:p) {
            for(j in 1:p) {
                  if(i != j) K_Mm[i, j] <- -1/p
            }
      }
      Sigma_Mm <- solve(K_Mm)

      X_Mm <- rmvnorm(n, rep(0,p), Sigma_Mm)
      MH_trunc_X_Mm <- MH(K_Mm, X_Mm, nMH, burn_rate = 0.1, plot_dens = FALSE, statistics = FALSE)
      time_p_fin <- Sys.time()
      time_p_vector_MH <- c(time_p_vector_MH, difftime(time_p_fin, time_p_ini, units = "mins"))
}

data_time <- as_data_frame(cbind(p_vector, time_p_vector_gibbs, time_p_vector_MH))
colnames(data_time) <- c("p", "3.1", "3.2")
data_time_melt <- melt(data_time, id.vars = "p")
colnames(data_time_melt) <- c("p", "Algorithm", "value")
plot_time <- ggplot(data_time_melt, aes(x = p, y = value, group = Algorithm, colour = Algorithm)) +
      geom_line(lwd = 2) + geom_point() + ylab("minutes") + theme_bw()


# Data sets exploration
data("math")
X_mardia <- as.matrix(math)
K_mardia <- solve(cor(X_mardia))

MH_X_mardia <- MH(K_mardia, X_mardia, nMH, burn_rate = 0.1)
MH_X_mardia$dens_plot_1
MH_X_mardia$dens_plot_2
MH_X_mardia$Gibbs_summary$mean
MH_X_mardia$is_m_matrix_prop

X_mardia <- as.matrix(math)[, -4]
K_mardia <- solve(cor(X_mardia))

MH_X_mardia <- gibbs_trunc(K_mardia, X_mardia, nMH, burn_rate = 0.1)
MH_X_mardia$dens_plot_1
MH_X_mardia$dens_plot_2
MH_X_mardia$Gibbs_summary$mean
MH_X_mardia$is_m_matrix_prop

MH_X_mardia <- MH(K_mardia, X_mardia, nMH, burn_rate = 0.1)
MH_X_mardia$dens_plot_1
MH_X_mardia$dens_plot_2
MH_X_mardia$Gibbs_summary$mean
MH_X_mardia$is_m_matrix_prop

latex_mardia_corr <- xtable(solve(cor(X_mardia)), 
                       align=rep("", ncol(cor(X_mardia)) + 1))
print(latex_mardia_corr, floating=FALSE,
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

latex_mardia_summary <- xtable(MH_X_mardia$Gibbs_summary$mean, 
                            align = rep("", ncol(cor(X_mardia)) + 1),
                            digits = 4)

print(latex_mardia_summary, floating=FALSE,
      hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)
