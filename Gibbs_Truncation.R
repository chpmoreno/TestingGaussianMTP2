gibbs_trunc <- function(K_ini, X, ngibbs, Gamma_sign = NULL, burn_rate = 0.1, plot_dens = TRUE,
                        statistics = TRUE) {
      S <- t(X) %*% X
      n <- nrow(X)
      p <- ncol(X)
      burn   <- ngibbs * burn_rate
      if(is.null(Gamma_sign)) {
            Gamma_sign = sign(K_ini)
      }
      K <- array(NA, c(p, p, ngibbs + 1))
      K[ , , 1] <- K_ini
      n_hold <- 0
      for(b in 2:(ngibbs + 1)) {
            K[ , , b] <- K[ , , b - 1]
            for(i in 1:(p-1)) {
                  for(j in (i+1):p) {
                        if(Gamma_sign[i, j] == 0) {
                              K[i, i, b] <- rWishart(1, n + p + 1, solve(as.matrix(S[i, i])))
                              K[j, j, b] <- rWishart(1, n + p + 1, solve(as.matrix(S[j, j])))
                              K[i, j, b] <- 0
                              K[j, i, b] <- 0
                        } else{
                              R <- setdiff(1:p, c(i, j))
                              rhs <- -K[c(i, j), R, b] %*% solve(K[R, R, b]) %*% K[R, c(i, j), b]
                              repeat {
                                    n_hold <- n_hold + 1
                                    tK <- rWishart(1, n + p + 1, solve(diag(2) + S[c(i, j), c(i, j)]))[ , , 1] 
                                    if (Gamma_sign[i, j] * tK[1, 2] >= Gamma_sign[i, j] * rhs[1, 2]) break;
                              }
                              K[c(i, j), c(i, j), b] <- tK - rhs
                        }
                  }
            }
      }
      K <- K[ , , -(1:(burn + 1))]
      
      if(plot_dens) {
            K_df_names <- NULL
            K_df <- NULL
            for(i in 1:p) {
                  for(j in 1:p) {
                        if(i <= j){
                              if(sum(K[i, j, ]) != 0) {
                                    K_df_names = c(K_df_names, paste0("K[", i, ",", j, "]"))
                                    K_df <- cbind(K_df, K[i, j, ])
                              }
                        }
                  }
            }
            K_df <- as_data_frame(K_df)
            colnames(K_df) <- K_df_names
            K_df_melt <- melt(K_df)
            
            if(ncol(K_df) < 6) {
                  dens_plot_1 <- ggplot(data = K_df_melt, aes(value, colour = variable, fill = variable, group = variable)) +
                        geom_density(alpha = 0.9) + geom_hline(yintercept=0, colour="gray40", size = 0.8) + 
                        facet_grid(~variable, scales = "free_x") + theme_bw()
            } else {
                  dens_plot_1 <- ggplot(data = K_df_melt, aes(value, colour = variable, fill = variable, group = variable)) +
                        geom_density(alpha = 0.9) +  geom_hline(yintercept=0, colour="gray40", size = 0.8) + 
                        facet_wrap(~variable, ncol = ncol(K_df) / 5) + theme_bw()
            }
            
            dens_plot_2 <- ggplot(data = K_df_melt, aes(value, colour = variable, fill = variable, group = variable)) +
                  geom_density(alpha = 0.9) +  geom_hline(yintercept=0, colour="gray40", size = 0.8) +
                  theme_bw()      
      } else {
            dens_plot_1 <- NULL
            dens_plot_2 <- NULL
      }
      
      if(statistics) {
            Gibbs_summary <- list()
            Gibbs_summary[[1]] <- apply(K, 1:2, mean)
            Gibbs_summary[[2]] <- apply(K, 1:2, sd)
            Gibbs_summary[[3]] <- apply(K, 1:2, min)
            Gibbs_summary[[4]] <- apply(K, 1:2, max)
            names(Gibbs_summary) <- c("mean", "sd", "min", "max")
            
            is_m_matrix_vector = NULL
            for(i in 1:dim(K)[3]) {
                  is_m_matrix_vector <- c(is_m_matrix_vector, is_m_matrix(K[ , , i]))
            }
            # Proportion of M-Matrices in the Gibbs Sampler
            is_m_matrix_prop = sum(is_m_matrix_vector) / dim(K)[3]
      } else {
            Gibbs_summary <- NULL
            is_m_matrix_prop <- NULL
      }
      
      return(list(K = K, Gibbs_summary = Gibbs_summary, 
                  dens_plot_1 = dens_plot_1, dens_plot_2 = dens_plot_2,
                  n_hold = n_hold, is_m_matrix_prop = is_m_matrix_prop))
}