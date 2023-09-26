#' Align the estimated clusters and the ground truth clusters
#'
#' @param est A numeric vector; the clusters.
#' @param truth A numeric vector; the clusters.
#' @param force_merge TRUE or FALSE; whether the merge the clusters in `est` when
#' the number of clusters does not match the ground truth.
#'
#' @return
#' @keywords internal
align_clusters <- function(est, truth, force_merge = FALSE) {
  tbl <- table(est, truth)
  est_labels <- as.numeric(rownames(tbl))
  truth_labels <- as.numeric(colnames(tbl))

  mapped_est <- est
  if (force_merge) {
    mapped_clusters <- apply(as.matrix(tbl), 1, which.max)
    for (i in seq_along(mapped_clusters)) {
      ind <- mapped_clusters[i]
      mapped_est[est == est_labels[i]] <- truth_labels[ind]
    }
  } else {
    rev_mapped_clusters <- apply(as.matrix(tbl), 2, which.max)
    for (i in seq_along(rev_mapped_clusters)) {
      ind <- rev_mapped_clusters[i]
      mapped_est[est == est_labels[ind]] <- truth_labels[i]
    }
    mapped_est[!(est %in% est_labels[rev_mapped_clusters])] <- -1
  }
  mapped_est
}


#' Compare the estimated (and aligned) clusters with the ground truth clusters
#'
#' @param est A numeric vector; the clusters.
#' @param truth A numeric vector; the clusters.
#' @param force_merge TRUE or FALSE; whether the merge the clusters in `est` when
#' the number of clusters does not match the ground truth.
#'
#' @return
#' @keywords internal
compare_clusters <- function(est, truth, force_merge = FALSE) {
  sum(align_clusters(est, truth, force_merge) == truth) * 100 / length(truth)
}


# Testing
testthat::test_that("Test align_clusters", {
  testthat::expect_true(all(
    align_clusters(c(1, 1, 1, 2, 2, 3, 7, 8, 9),
                   c(4, 4, 4, 5, 5, 6, 4, 4, 4),
                   TRUE) == c(4, 4, 4, 5, 5, 6, 4, 4, 4)
  ))

  testthat::expect_true(all(
    align_clusters(c(3, 3, 4, 4, 5, 2, 3, 4),
                   c(2, 2, 3, 3, 4, 2, 2, 2),
                   TRUE) == c(2, 2, 3, 3, 4, 2, 2, 3)
  ))

  testthat::expect_true(all(
    align_clusters(c(1, 1, 1, 2, 2, 3, 7, 8, 9),
                   c(4, 4, 4, 5, 5, 6, 4, 4, 4),
                   FALSE) == c(4, 4, 4, 5, 5, 6, -1, -1, -1)
  ))

  testthat::expect_true(all(
    align_clusters(c(3, 3, 4, 4, 5, 2, 3, 4),
                   c(2, 2, 3, 3, 4, 2, 2, 2),
                   FALSE) == c(2, 2, 3, 3, 4, -1, 2, 3)
  ))
})


testthat::test_that("Test compare_clusters", {
  testthat::expect_equal(
    compare_clusters(c(1, 1, 1, 2, 2, 3, 7, 8, 9),
                     c(4, 4, 4, 5, 5, 6, 4, 4, 4),
                     TRUE),
    100  # Perfect
  )

  testthat::expect_equal(
    compare_clusters(
      c(1, 1, 1, 2, 2, 3, 7, 8, 9),
      c(4, 4, 4, 5, 5, 6, 4, 4, 4),
      FALSE),
    100 * 6 / 9  # Misses 3
  )
})
