set.seed(1122)
N = 300
U = 200

hc_acc <- hc_ARI <- hc_NMI <- NULL
K = 5
for(RDR_var in seq(0.1, 2, 0.2)){
    sims <- get_sim_data(K = K, N = N, U = U, expr = F, RDR_var = RDR_var, RDR_outlier_cnt = 0, CNV_overlap= FALSE)
    
    hc <- hclust(dist(t(sims$RDR)), method = 'ward.D2')
    K_hat <- K + 5
    I_label_hc <- cutree(hc, k = K_hat)
    fin_hc <- align_clusters(I_label_hc, sims$cluster_true)
    hc_acc <- c(hc_acc, sum(fin_hc == sims$cluster_true)*100/N)
    hc_ARI <- c(hc_ARI, aricode::ARI(I_label_hc, sims$cluster_true)*100)
    hc_NMI <- c(hc_NMI, aricode::NMI(I_label_hc, sims$cluster_true)*100)
}
## not directly comparable because the simulated data is different but the trend is going up 
## with increasing K, for all criterion
print(hc_acc)
print(hc_ARI)
print(hc_NMI)