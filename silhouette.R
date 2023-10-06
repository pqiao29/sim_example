#### silhouette score 
get_silhouette <- function(clust_est){
    silhouette <- rep(NA, N)
    for(i in 1:N){
        tmp_k <- clust_est[i]
        same_cluster <- clust_est == tmp_k
        
        a <- ifelse(sum(same_cluster) > 1, sum(profile_dist[i, same_cluster])/(sum(same_cluster) - 1), 0)
        b <- NULL
        for(k in (1:length(unique(clust_est)))[-tmp_k]){
            b <- c(b, mean(profile_dist[i, clust_est == k]))
        }
        b <- min(b)
        
        if(a > 0 && b == 0) silhouette[i] <- -1
        if(a == 0 && b == 0) silhouette[i] <- 0
        if(b > 0) silhouette[i] <- (b - a)/max(a, b)
    }
    silhouette
}


## should decrease
hc_silhouette <- NULL
RDR_var <- 1
for(K in 3:6){
    sims <- get_sim_data(K = K, N = N, U = U, expr = F, RDR_var = RDR_var, RDR_outlier_cnt = 0, CNV_overlap= FALSE)
    
    #### Build distance
    profile_dist <- matrix(NA, N, N) 
    diag(profile_dist) <- 0
    for(i1 in 1:(N - 1)){
        for(i2 in (i1 + 1):N){
            profile_dist[i1, i2] <- profile_dist[i2, i1] <- sum(sims$cell_level_states[, i1] != sims$cell_level_states[, i2])
        }
    }
    
    #### hclust
    hc <- hclust(dist(t(sims$RDR)), method = 'ward.D2')
    I_label_hc <- cutree(hc, k = K_hat)
    
    hc_silhouette <- c(hc_silhouette, mean(get_silhouette(I_label_hc)))
    
}
print(hc_silhouette)

## should decrease
K <- 5
for(RDR_var in seq(0.1, 2, 0.2)){
    sims <- get_sim_data(K = K, N = N, U = U, expr = F, RDR_var = RDR_var, RDR_outlier_cnt = 0, CNV_overlap= FALSE)
    
    #### Build distance
    profile_dist <- matrix(NA, N, N) 
    diag(profile_dist) <- 0
    for(i1 in 1:(N - 1)){
        for(i2 in (i1 + 1):N){
            profile_dist[i1, i2] <- profile_dist[i2, i1] <- sum(sims$cell_level_states[, i1] != sims$cell_level_states[, i2])
        }
    }
    
    #### hclust
    hc <- hclust(dist(t(sims$RDR)), method = 'ward.D2')
    I_label_hc <- cutree(hc, k = K_hat)
    
    hc_silhouette <- c(hc_silhouette, mean(get_silhouette(I_label_hc)))
    
}
print(hc_silhouette)
