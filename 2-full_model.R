## Fig 2C: state accuracy (N = 200, U = 300, iterations: 500+500)

K = 5                 ## True cluster number
N = 100               ## cell counts
U = 200               ## gene counts
RDR_var = 1           ## variance of RDR (x-axis of Fig 2C)
CNV_overlap = F       ## If true, CN profiles between different clones will be more similar, harder to identify K and clustering 
Q_to_neutral = 4      ## In generating CN profiles: for each gene to the next, the probability of jumping from a non-neutral to a neutral state is "Q_to_neutral" times larger than jumping to another non-neutral state

### simulate RDR
sims <- get_sim_data(K = K, N = N, U = U, expr = F, RDR_var = RDR_var, 
                     Q_to_neutral = Q_to_neutral, CNV_overlap = CNV_overlap)

plot_inout(sims$RDR, "RDR", list(sims$cluster_true), CN_states = sims$states_true)

### fit model 
res <- Chloris(sims$RDR, K = 10, 
               Gibbs_tol = 100, burnin_tol = 100,
               cluster_shrink_tol = 20, min_cluster_size = 2)
plot_inout(sims$RDR, "RDR", cluster_labels = list(res$cluster_est, sims$cluster_true),
           CN_states = res$state_est, state_mean = apply(res$par_record[, , 1], 2, mean))


### Align resulting labels with true labels to enable comparison
import::from("~/Documents/Projects/old_simulations/res_alignment.R", cluster_align2, state_align2, align_fin3)
fin_RDR <- align_fin3(est_I_label = res$cluster_est, est_states = res$state_est,
                      true_I_label = sims$cluster_true, verbose = FALSE)

cat("\nClustering accuracy: ", sum(fin_RDR$clusters_est == sims$cluster_true)*100/length(fin_RDR$clusters_est), "%\n")
cat("State accuary:", sum(fin_RDR$states_est == sims$states_true[sort(unique(fin_RDR$cluster_alignment)), ])*100/length(sims$states_true[sort(unique(fin_RDR$cluster_alignment)), ]), "%")

