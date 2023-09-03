### Supplementary Fig 2

S = 3        
N = 100 
U = 200

sims <- get_sim_data(RDR = F, BAF = T, K = 6, N = N, U = U, S = S,
                     BAF_missing_percent = 0.9)


### plot truth
plot_inout(sims$A/sims$D, "BAF", list(sims$cluster_true),
           CN_states = sims$states_true)


res <- Chloris(A = sims$A, D = sims$D, S = S, init = "random",
               Gibbs_tol = 100, burnin_tol = 100)

### plot result
plot_inout((sims$A/sims$D), "BAF", list(res$cluster_est, sims$cluster_true), 
           CN_states = res$state_est)

### get accuracy 
norm_clone <- which.min(apply(res$state_est, 1, function(x) sum(x != S)))
group1 <- sims$cluster_true == 1
group2 <- res$cluster_est == norm_clone
plot_inout((sims$A/sims$D)[keep_snps, ], "BAF", list(group2, group1))

cat("Total accuracy:", sum(group1 == group2)*100/N, "%")
cat("Purity: ", sum(group1 & group2)*100/sum(group2), "%")
