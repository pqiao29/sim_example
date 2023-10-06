cluster <- function(n = 100, id, mu, sigma) {
  if (missing(id)) id <- sample(1e3, 1)
  if (missing(mu)) mu <- rnorm(2, 0, 4)
  if (missing(sigma)) sigma <- rgamma(2, 4, 3)

  x <- rnorm(n, mu[1], sigma[1])
  y <- rnorm(n, mu[2], sigma[2])
  data.frame(x = x, y = y, id = id)
}


set.seed(20231006)
seeds <- sample(1e8, 1e3)
set.seed(seeds[1])
true_K <- 4
df0 <- do.call(rbind, lapply(1:true_K, \(id) cluster(100, id)))
plot(df0$x, df0$y, col = df0$id, pch = 19)


for (k in 5:10) {
  hcl <- df0[c("x", "y")] |>
    dist() |>
    hclust(method = "ward.D2") |>
    cutree(k = k)
  hcl

  perf <- compare_clusters(hcl, df0$id, force_merge = FALSE)
  message(sprintf("k=%s; accuracy: %s", k, perf))
}
