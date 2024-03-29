state_align <- function(est, true, par_record, RDR_outlier_cnt){
  
  S = 4
  RDR_levels <- log2(c(0.5, 1, 1.5, 2))
  state_est <- est
  for(s in 1:S){
    state_label <- which.min(abs(RDR_levels - mean(par_record[, s, 1])))
    if(state_label != s) state_est[est == s] <- state_label
  } 
  
  
  K_est <- nrow(state_est)
  K_true <- nrow(true)
  state_dist <- matrix(0, K_est, K_true)
  for(k_est in 1:K_est){
    for(k_true in 1:K_true){
      state_dist[k_est, k_true] <- sum(true[k_true, ] != state_est[k_est, ])
    }
  }
    
    100 - sum(apply(state_dist, 2, min))*100/(length(state_est))
}