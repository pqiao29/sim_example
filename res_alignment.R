rm_empty_cluster <- function(est){
    
    empty_clusters <- which(colSums(est$I) == 0)
    if(length(empty_clusters) > 0){
        est$H <- est$H[-empty_clusters, ]
        est$I <- est$I[, -empty_clusters]
        est$I_label <- apply(est$I, 1, which.max)
        est$K <- ncol(est$I)
        message(paste0("Removing ", length(empty_clusters), " empty clusters ... \n"))
    }
    
    return(est)
}

label_align_heuristic <- function(rematch){
    
    S_est <- nrow(rematch)
    S_true <- ncol(rematch)
    
    state_reorder <- rep(NA, S_est)
    unmatched_states_est <- 1:S_est
    unmatched_states_true <- 1:S_true
    while(length(unmatched_states_est) > 0 && length(unmatched_states_true) > 0){
        idx <- which(rematch[unmatched_states_est, unmatched_states_true, drop = FALSE] == max(rematch[unmatched_states_est, unmatched_states_true, drop = FALSE]), arr.ind = TRUE)
        if(nrow(idx) > 1){
            if(length(unique(idx[, 1])) > 1 && length(unique(idx[, 2])) > 1){ ## maximum entries appear on different rows and columns
                tmp_flag <- 1
                found <- FALSE
                while(!found && tmp_flag <= nrow(idx)){
                    if(sum(idx[, 1] == idx[tmp_flag, 1]) == 1 && sum(idx[, 2] == idx[tmp_flag, 2]) == 1){
                        found <- TRUE  
                        idx <- idx[tmp_flag, ]
                    }else{
                        tmp_flag<- tmp_flag + 1
                    }
                }
            }else{ ### maximum entries either on the same row or same column
                diff_dim <- ifelse(length(unique(idx[, 1])) == 1, "col", "row")
                tmp_rematch <- rematch[unmatched_states_est, unmatched_states_true, drop = FALSE]
                tmp_idx <- rep(NA, 2)
                if(diff_dim == "col"){
                    tmp_idx[1] <- unique(idx[, 1])
                    tmp_rematch <- tmp_rematch[-unique(idx[, 1]), idx[, 2], drop = FALSE]
                    tmp_idx[2] <- which.min(colSums(tmp_rematch))
                }else{
                    tmp_idx[2] <- unique(idx[, 2])
                    tmp_rematch <- tmp_rematch[idx[, 1], -unique(idx[, 2]), drop = FALSE]
                    tmp_idx[1] <- which.min(rowSums(tmp_rematch))
                }
                idx <- tmp_idx
            }
        }
        state_reorder[unmatched_states_est[idx[1]]] <- unmatched_states_true[idx[2]]
        unmatched_states_est <- unmatched_states_est[-idx[1]]
        unmatched_states_true <- unmatched_states_true[-idx[2]]
    }
    
    if(S_est > S_true){
        for(s in unmatched_states_est){
            state_reorder[s] <- which.max(rematch[s, ])
        }
    }
    
    if(S_est <= S_true){
        state_reorder[is.na(state_reorder)] <- unmatched_states_true[which.max(rematch[unmatched_states_est, unmatched_states_true])] 
    }else{ 
        if(length(unmatched_states_true) == 1){
            idx <- unmatched_states_est[which.max(rematch[unmatched_states_est, unmatched_states_true])]
            state_reorder[idx] <- unmatched_states_true  ## !!! !!! possible bug when S_est - S_true >=2 !!! !!!  
            state_reorder[is.na(state_reorder)] <- S_est
        }
    }
    
    return(state_reorder)
}  ## heuristic algorithm to find the permutation that maximizes the trace

label_align_exhaust <- function(est, true, S1, S2, verbose = TRUE){
    #### create all distinct permutations
    p <- list()
    p[1:S1] <- rep(list(1:S2), S1)
    all <- expand.grid(p, stringsAsFactors = FALSE) 
    perms <- as.matrix(all[apply(all, 1, function(x) {length(unique(x)) == S1}),])
    dimnames(perms) <- NULL
    #### Compute similarity for each permutation
    similarity <- rep(NA, nrow(perms))
    if(is.vector(est)){
        N <- length(true)
        for(perm_idx in 1:nrow(perms)){  
            tmp_aligned <- as.integer(as.character(factor(est, levels = 1:S1, labels = perms[perm_idx, ])))
            similarity[perm_idx] <- sum(tmp_aligned == true)
        }
    }else{
        if(is.matrix(est)){
            N <- prod(dim(true))
            for(perm_idx in 1:nrow(perms)){  
                tmp_aligned <- t(apply(est, 1, function(x) as.integer(as.character(factor(x, levels = 1:S1, labels = perms[perm_idx, ])))))
                similarity[perm_idx] <- sum(tmp_aligned == true)
            }
        }else{
            stop("The est argument in label_align_exhaust() needs to be either a vector or a matrix!")
        }
    }
    #### Find permutation that maximizes similarity
    reorder_idx <- which(similarity == max(similarity))
    identical_idx <- which(apply(perms, 1, function(x) identical(x, 1:S1)))
    #### return only if label has been switched
    if(!(identical_idx %in% reorder_idx)){
        if(verbose) cat("\nSimilarity with flag before switching: ", similarity[identical_idx]*100/N, 
                        "% \nSimilarity with flag after switching: ",  max(similarity)*100/N, "%.\n")
        if(length(reorder_idx) > 1 && verbose) cat(length(reorder_idx), " best match detected!!! Taking the first one. \n")
        reorder <- perms[reorder_idx[1], ]
        return(reorder)
    }else{
        return(NULL)
    }
}

cluster_align2 <- function(est, true_I, fin = FALSE, exhaust = FALSE, verbose){  ## true_I is similar to est$I_label
    
    if(fin){
        est_I <- est
    }else{
        est_I <- est$I_label
    }
    
    ### Prepare
    K1 <- length(unique(est_I))
    K2 <- length(unique(true_I))
    if(max(est_I) > K1) est_I <- as.integer(factor(est_I, levels = unique(est_I), labels = 1:K1))  ## in case some clusters are empty, e.g. 2, 3, 4, 3, 3 -> 1, 2, 3, 2, 2
    if(K1 > K2 && exhaust){
        exhaust <- FALSE
        if(verbose) cat("Estimates got more clusters than true/flag, going for heuristic algorithm\n")
    }
    
    ### Find optimal order
    if(!exhaust){
        rematch <- matrix(NA, K1, K2)
        for(k1 in 1:K1){
            for(k2 in 1:K2){
                rematch[k1, k2] <-  sum((est_I == k1)*(true_I == k2)) 
            }
        } 
        cluster_reorder <- label_align_heuristic(rematch)
        if(identical(cluster_reorder, 1:K1)){
            cluster_reorder <- NULL
        }else{
            tmp_aligned_I <- as.integer(as.character(factor(est_I, levels = 1:K1, labels = cluster_reorder)))
            similarity_after <- sum(tmp_aligned_I == true_I)
            similarity_before <- sum(est_I == true_I)
            if(similarity_after <= similarity_before){
                cluster_reorder <- NULL
            }else{
                if(verbose)  cat("Cluster label switching (heuristic)! \nSimilarity with flag before switching: ", similarity_before*100/length(true_I), 
                                 "% \nSimilarity with flag after switching: ",  similarity_after*100/length(true_I), "%.\n")
            }
        }
    }else{
        if(verbose)  cat("\nAligning clusters: \n")
        cluster_reorder <- label_align_exhaust(est_I, true_I, K1, K2, verbose = verbose)
    }
    
    ### Wrap up
    if(!is.null(cluster_reorder)){
        aligned_I <- as.integer(as.character(factor(est_I, levels = 1:K1, labels = cluster_reorder)))
        if(!fin){
            
            reverse_order <- rep(NA, est$K)
            for(target_idx in cluster_reorder){
                current_idx <- est_I[aligned_I == target_idx][1]
                reverse_order[target_idx] <- current_idx
            }
            if(any(is.na(reverse_order))) reverse_order[is.na(reverse_order)] <- (1:est$K)[-reverse_order[!is.na(reverse_order)]]
            
            est$I <- est$I[, reverse_order]
            est$Q <- est$Q[reverse_order]
            est$H <- est$H[reverse_order, ]
            est$I_label <- aligned_I
            est$cluster_switch_signal <- TRUE
            return(est) 
        }else{
            return(list("order" = cluster_reorder, "cluster" = aligned_I, "switched" = TRUE)) 
        }
    }else{
        if(!fin){
            est$cluster_switch_signal <- FALSE
            return(est)
        }else{
            return(list("order" = 1:K1, "cluster" = est_I, "switched" = FALSE))
        }
    }
}

state_align2 <- function(est, true_H, fin = FALSE, exhaust = FALSE, verbose){
    
    if(fin){
        if(!("matrix" %in% class(est))) stop("If fin == TRUE, est needs to be a matrix!")
        est_H <- est
        S1 <- max(est_H)
    }else{
        if(!("list" %in% class(est))) stop("if fin == FALSE, est needs to be a list!")
        est_H <- est$H
        S1 <- est$S
    }
    
    S2 <- max(true_H)
    K <- nrow(true_H)
    if(K != nrow(est_H)) stop("In state alignment, the numbers of clusters of est and true don't agree")
    if(S1 > S2 && exhaust){
        exhaust <- FALSE
        if(verbose) cat("Estimates got more states than true/flag, going for heuristic algorithm\n")
    }
    
    if(!exhaust){
        state_rematch <- matrix(0, S1, S2)
        for(s1 in 1:S1){
            for(s2 in 1:S2){
                for(k in 1:K){
                    state_rematch[s1, s2] <- state_rematch[s1, s2] + sum((est_H[k, ] == s1)*(true_H[k, ] == s2))   
                }
                
            }
        } 
        state_relabel <- label_align_heuristic(state_rematch)
        if(identical(state_relabel, 1:S1)){   ### check if state_relabel doesn't improve alignment, then state_relabel <- NULL
            state_relabel <- NULL
        }else{
            tmp_aligned_H <- t(apply(est_H, 1, function(x) as.integer(as.character(factor(x, levels = 1:S1, labels = state_relabel)))))
            similarity_after <- sum(tmp_aligned_H == true_H)
            similarity_before <- sum(est_H == true_H)
            if(similarity_after <= similarity_before){
                state_relabel <- NULL
            }else{
                if(verbose)  cat("State label switching (heuristic)! \nSimilarity with flag before switching: ", similarity_before*100/prod(dim(true_H)), 
                                 "% \nSimilarity with flag after switching: ",  similarity_after*100/prod(dim(true_H)), "%.\n")
            }
        }
    }else{
        if(verbose)  cat("\nAligning states: \n")
        state_relabel <- label_align_exhaust(est_H, true_H, S1, S2, verbose = verbose)
    }
    
    if(!is.null(state_relabel)){
        aligned_H <- t(apply(est_H, 1, function(x) as.integer(as.character(factor(x, levels = 1:S1, labels = state_relabel)))))
        if(!fin){
            est$mu <- est$mu[state_relabel]
            est$sigma <- est$sigma[state_relabel]
            est$theta <- est$theta[state_relabel]
            est$Q <- lapply(est$Q, function(x) x[state_relabel, state_relabel])
            est$pi <- lapply(est$pi, function(x) x[, state_relabel])
            est$H <- aligned_H
            est$state_switch_signal <- TRUE 
            return(est)
        }else{
            return(list("order" = state_relabel, "state" = aligned_H, "switched" = TRUE))
        }
    }else{
        if(!fin){
            est$state_switch_signal <- FALSE
            return(est)
        }else{
            return(list("order" = 1:S1, "state" = est, "switched" = FALSE))
        }
    }
    
}


align_fin2 <- function(est, I_label_true, states_true, verbose = TRUE){

    est <- rm_empty_cluster(est)

    ##### Align clusters
    cluster_order <- cluster_align2(est$I_label, I_label_true, fin = TRUE, exhaust = TRUE, verbose = verbose)
    
    ##### Align states
    H_est <- est$H[sort(unique(est$I_label)), ] ## get rid of CNV profiles of empty clusters
    S1 <- length(unique(c(H_est)))
    if(max(H_est) > S1){
        dim_H_est <- dim(H_est)
        H_est <- as.numeric(as.character(factor(H_est, levels = unique(c(H_est)), labels = 1:S1))) ## get rid of empty states
        H_est <- matrix(H_est, dim_H_est[1], dim_H_est[2])
    } 
    state_order <- state_align2(H_est, states_true[cluster_order$order, ], fin = TRUE, exhaust = TRUE, verbose = verbose)
    
    res <- list()
    res$clusters_est <- cluster_order$cluster
    res$clusters_true <- I_label_true
    res$cluster_alignment <- cluster_order$order
    res$states_est <- state_order$state

    return(res)
}


### No label switching considered, only to align with truth for assessment
align_fin3 <- function(est_I_label, est_states, true_I_label, verbose = TRUE){
    
    ##### Align clusters
    cluster_order <- cluster_align2(est_I_label, true_I_label, fin = TRUE, exhaust = TRUE, verbose = verbose)
    
    ##### Align states (assume no cluster label switch during sampling)
    non_empty_clusters <- sort(unique(est_I_label))
    H_est <- est_states[non_empty_clusters, ]
    identified_Ks <- unique(cluster_order$order)
    
    H_est_aligned <- matrix(NA, length(identified_Ks), ncol(H_est))
    if(length(unique(est_I_label)) >= length(unique(true_I_label))){
        for(k in identified_Ks){
            tmp_H_k <- H_est[cluster_order$order == k, , drop = FALSE]
            if(sum(cluster_order$order == k) > 1){
                states_concensus_k <- apply(tmp_H_k, 2, function(x) {
                    tmp <- table(x)
                    tmp_idx <- which(tmp == max(tmp))
                    if(length(tmp_idx) == 1) ret <- tmp_idx else ret <- sample(tmp_idx, 1) 
                    as.integer(names(ret))
                })
            }else states_concensus_k <- tmp_H_k
            H_est_aligned[k, ] <- states_concensus_k
        }
    }else{
        identified_Ks <- sort(identified_Ks)
        H_idx <- 1
        for(k in identified_Ks){
            tmp_H_k <- H_est[cluster_order$order == k, , drop = FALSE]
            if(sum(cluster_order$order == k) > 1){
                states_concensus_k <- apply(tmp_H_k, 2, function(x) {
                    tmp <- table(x)
                    tmp_idx <- which(tmp == max(tmp))
                    if(length(tmp_idx) == 1) ret <- tmp_idx else ret <- sample(tmp_idx, 1) 
                    as.integer(names(ret))
                })
            }else states_concensus_k <- tmp_H_k
            H_est_aligned[H_idx, ] <- states_concensus_k
            H_idx <- H_idx + 1
        }
    }
    
       
    res <- list()
    res$clusters_est <- cluster_order$cluster
    res$clusters_true <- true_I_label
    res$cluster_alignment <- cluster_order$order
    res$states_est <- H_est_aligned
    
    return(res)
}
