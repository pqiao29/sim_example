## Fig 2B Clustering accuracy, compare with hclust

K = 5               ## True cluster number
K_hat = K + 1       ## Input cluster number, same for Chloris and hclust
RDR_outlier_cnt = 5 ## outlier cell counts. Outliers cell: random CN profiles, 4 times variance than RDR_var
N = 100             ## (non-outlier) cell counts
U = 200             ## gene counts
RDR_var = 1       ## RDR variance

### Simulate RDR
sims <- get_sim_data(K = K, N = N, U = U, expr = F, RDR_var = RDR_var, RDR_outlier_cnt = RDR_outlier_cnt)
plot_inout(sims$RDR, "RDR", list(c(sims$cluster_true, rep(K + 1, RDR_outlier_cnt))), CN_states = rbind(sims$states_true, rep(2, U)), lim = c(-3, 3))

### Run model
res <- Chloris(sims$RDR, burnin_tol = 100, Gibbs_tol = 100, cluster_shrink_tol = NULL, K = K_hat, init = "random")
res_RDR_aligned <- res

### Align resulting cluster labels with true cluster labels, to enable comparison
import::from("~/Documents/Projects/old_simulations/res_alignment.R", cluster_align2, state_align2, align_fin3)
fin_RDR <- align_fin3(est_I_label = res_RDR_aligned$cluster_est, est_states = res_RDR_aligned$state_est, 
                      true_I_label = c(sims$cluster_true, rep(K + 1, RDR_outlier_cnt)), verbose = FALSE)

### Remove the cluster formed by outliers
if((K + 1) %in% unique(fin_RDR$cluster_alignment)){ ## outlier cluster identified
    fin_RDR$states_est <- fin_RDR$states_est[-nrow(fin_RDR$states_est), ]
    fin_RDR$cluster_alignment <- fin_RDR$cluster_alignment[fin_RDR$cluster_alignment != (K + 1)]
}

cat("\nCluster accuracy\n")
cat("from Chloris:", sum(fin_RDR$clusters_est[1:N] == sims$cluster_true[1:N])*100/N, "%\n")

### Run hclust
hc <- hclust(dist(t(sims$RDR)), method = 'ward.D2')
I_label_hc <- cutree(hc, k = K_hat)
fin_hc <- cluster_align2(I_label_hc[1:N], sims$cluster_true, fin = TRUE, exhaust = TRUE, verbose = FALSE)
cat("from hclust:", sum(fin_hc$cluster == sims$cluster_true)*100/length(sims$cluster_true), "%")

