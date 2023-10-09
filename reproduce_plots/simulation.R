run <- function(RDR_var, K, seed, burnin_tol, Gibbs_tol){
  
  set.seed(seed)
  ret <- matrix(NA, 3, 2)
  colnames(ret) <- c("Clusters", "States")
  rownames(ret) <- c("RDR", "both", "hclust")
  
  K_hat = round(K*2.5)       ## Input cluster number, same for Chloris and hclust
  N = 300
  U = 250
  
  sims <- get_sim_data(K = K, N = N, U = U, expr = F, RDR_var = RDR_var, RDR_outlier_cnt = 0, CNV_overlap= FALSE)
  
  #### Chloris
  res_RDR <- Chloris(sims$RDR, burnin_tol = burnin_tol, Gibbs_tol = Gibbs_tol, cluster_shrink_tol = 20, min_cluster_size = 1, K = K_hat, init = "hclust", verbose = F)
  res_both <- Chloris(sims$RDR, A = sims$A, D = sims$D, burnin_tol = burnin_tol, Gibbs_tol = Gibbs_tol, 
                      cluster_shrink_tol = 20, min_cluster_size = 1, K = K_hat, init = "hclust", verbose = F)
  
  #### hclust
  hc <- hclust(dist(t(sims$RDR)), method = 'ward.D2')
  I_label_hc <- cutree(hc, k = K_hat)
  
  
  ret[1, 2] <- state_align(res_RDR$state_est, sims$states_true, res_RDR$par_record, 0)
  ret[2, 2] <- state_align(res_both$state_est, sims$states_true, res_both$par_record, 0)

  ret[1, 1] <- sum(apply(table(res_RDR$cluster_est, sims$cluster_true), 2, max))*100/N
  ret[2, 1] <- sum(apply(table(res_both$cluster_est, sims$cluster_true), 2, max))*100/N
  ret[3, 1] <- sum(apply(table(I_label_hc, sims$cluster_true), 2, max))*100/N
  
  
  cat("\n")
  print(ret)
  ret
  
}