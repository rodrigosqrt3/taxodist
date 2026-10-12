# Clustering, ordination, heatmaps, and plot methods

test_that("plot.taxodist_cluster runs without error", {
  m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
              nrow = 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cl <- taxo_cluster(stats::as.dist(m))
  pdf(file = NULL)
  on.exit(dev.off())
  expect_no_error(plot(cl))
})

test_that("plot.taxodist_ord runs without error", {
  m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
              nrow = 3,
              dimnames = list(c("A","B","C"), c("A","B","C")))
  ord <- taxo_ordinate(stats::as.dist(m))
  pdf(file = NULL)
  on.exit(dev.off())
  expect_no_error(plot(ord))
  expect_invisible(plot(ord))
})

test_that("taxo_cluster returns correct S3 class", {
  mockery::stub(taxo_cluster, "distance_matrix", function(...) {
    m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
                nrow = 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
    stats::as.dist(m)
  })
  result <- taxo_cluster(c("A", "B", "C"), progress = FALSE)
  expect_s3_class(result, "taxodist_cluster")
})

test_that("taxo_cluster result contains hclust and dist", {
  mockery::stub(taxo_cluster, "distance_matrix", function(...) {
    m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
                nrow = 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
    stats::as.dist(m)
  })
  result <- taxo_cluster(c("A", "B", "C"), progress = FALSE)
  expect_s3_class(result$hclust, "hclust")
  expect_s3_class(result$dist, "dist")
})

test_that("taxo_cluster accepts a dist object directly", {
  m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
              nrow = 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  d <- stats::as.dist(m)
  result <- taxo_cluster(d)
  expect_s3_class(result, "taxodist_cluster")
})

test_that("taxo_cluster safely skips fewer than two taxa", {
  d <- stats::as.dist(matrix(
    0, nrow = 1, ncol = 1,
    dimnames = list("A", "A")
  ))
  expect_warning(
    result <- taxo_cluster(d),
    "At least two taxa"
  )
  expect_s3_class(result, "taxodist_cluster")
  expect_null(result$hclust)
  expect_equal(result$dist, d)
})

test_that("taxo_ordinate returns correct S3 class", {
  mockery::stub(taxo_ordinate, "distance_matrix", function(...) {
    m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
                nrow = 3,
                dimnames = list(c("A","B","C"), c("A","B","C")))
    stats::as.dist(m)
  })
  result <- taxo_ordinate(c("A", "B", "C"), progress = FALSE)
  expect_s3_class(result, "taxodist_ord")
})

test_that("taxo_ordinate result contains points, dist and GOF", {
  mockery::stub(taxo_ordinate, "distance_matrix", function(...) {
    m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
                nrow = 3,
                dimnames = list(c("A","B","C"), c("A","B","C")))
    stats::as.dist(m)
  })
  result <- taxo_ordinate(c("A", "B", "C"), progress = FALSE)
  expect_true(!is.null(result$points))
  expect_s3_class(result$dist, "dist")
  expect_true(!is.null(result$GOF))
})

test_that("taxo_ordinate points matrix has correct dimensions", {
  mockery::stub(taxo_ordinate, "distance_matrix", function(...) {
    m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
                nrow = 3,
                dimnames = list(c("A","B","C"), c("A","B","C")))
    stats::as.dist(m)
  })
  result <- taxo_ordinate(c("A", "B", "C"), k = 2, progress = FALSE)
  expect_equal(ncol(result$points), 2)
  expect_equal(nrow(result$points), 3)
})

test_that("taxo_ordinate accepts a dist object directly", {
  m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
              nrow = 3,
              dimnames = list(c("A","B","C"), c("A","B","C")))
  d <- stats::as.dist(m)
  result <- taxo_ordinate(d, k = 2)
  expect_s3_class(result, "taxodist_ord")
})

test_that("taxo_ordinate safely skips fewer than two taxa", {
  d <- stats::as.dist(matrix(
    0, nrow = 1, ncol = 1,
    dimnames = list("A", "A")
  ))
  expect_warning(
    result <- taxo_ordinate(d),
    "At least two taxa"
  )
  expect_s3_class(result, "taxodist_ord")
  expect_null(result$points)
})

test_that("taxo_ordinate reduces k for two taxa", {
  m <- matrix(c(0, 0.5, 0.5, 0), nrow = 2,
              dimnames = list(c("A", "B"), c("A", "B")))
  expect_warning(
    result <- taxo_ordinate(stats::as.dist(m), k = 2),
    "k.*reduced"
  )
  expect_equal(dim(result$points), c(2L, 1L))
})

test_that("taxo_ordinate rejects invalid k clearly", {
  m <- matrix(c(0, 0.5, 0.5, 0), nrow = 2)
  expect_error(
    taxo_ordinate(stats::as.dist(m), k = 0),
    "positive integer"
  )
})

test_that("summary.taxodist_ord computes variance and handles missing eigenvalues", {
  mock_ord <- structure(list(
    points = matrix(1:4, ncol = 2, dimnames = list(c("A", "B"), NULL)),
    dist   = stats::dist(1:2),
    GOF    = c(0.95, 0.95),
    eig    = c(2.0, 1.0, -0.5)
  ), class = "taxodist_ord")
  out <- capture.output(res <- summary(mock_ord))
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 2)
  expect_equal(res$Axis, c("PC1", "PC2"))
  expect_equal(round(res$Variance_Pct[1]), 67)
  expect_equal(round(res$Variance_Pct[2]), 33)
  mock_bad <- mock_ord
  mock_bad$eig <- NULL
  expect_message(res_bad <- summary(mock_bad), "Eigenvalues not found")
  expect_null(res_bad)
})

test_that("taxo_heatmap plots correctly and returns dist invisibly", {
  mock_dist <- stats::dist(matrix(1:4, ncol = 2))
  attr(mock_dist, "Labels") <- c("TaxonA", "TaxonB")
  pdf(file = NULL)
  res <- taxo_heatmap(mock_dist)
  dev.off()
  expect_s3_class(res, "dist")
})

test_that("taxo_heatmap safely skips fewer than two taxa", {
  d <- stats::as.dist(matrix(
    0, nrow = 1, ncol = 1,
    dimnames = list("A", "A")
  ))
  expect_warning(
    result <- taxo_heatmap(d),
    "At least two taxa"
  )
  expect_equal(result, d)
})

test_that("plot.taxodist_ord handles one-dimensional points", {
  ord <- structure(list(
    points = matrix(c(-0.25, 0.25), ncol = 1,
                    dimnames = list(c("A", "B"), "PC1")),
    dist = stats::dist(1:2),
    GOF = c(1, 1),
    eig = c(0.5, 0)
  ), class = "taxodist_ord")
  pdf(file = NULL)
  on.exit(dev.off())
  expect_invisible(plot(ord))
})

test_that("taxo_cluster handles NA in distance matrix gracefully", {
  m <- matrix(c(0, NA, NA, 0), nrow = 2, dimnames = list(c("A", "B"), c("A", "B")))
  d <- stats::as.dist(m)
  expect_warning(res <- taxo_cluster(d), "Distance matrix contains NA values")
  expect_null(res$hclust)
  expect_s3_class(res, "taxodist_cluster")
})

test_that("taxo_cluster handles infinite distances gracefully", {
  m <- matrix(c(0, Inf, Inf, 0), nrow = 2,
              dimnames = list(c("A", "B"), c("A", "B")))
  d <- stats::as.dist(m)
  expect_warning(
    res <- taxo_cluster(d),
    "infinite values"
  )
  expect_null(res$hclust)
  expect_equal(res$dist, d)
  expect_s3_class(res, "taxodist_cluster")
})

test_that("taxo_ordinate handles NA in distance matrix gracefully", {
  m <- matrix(c(0, NA, NA, 0), nrow = 2, dimnames = list(c("A", "B"), c("A", "B")))
  d <- stats::as.dist(m)
  expect_warning(res <- taxo_ordinate(d), "Distance matrix contains NA values")
  expect_null(res$points)
  expect_s3_class(res, "taxodist_ord")
})

test_that("taxo_ordinate handles infinite distances gracefully", {
  m <- matrix(c(0, Inf, Inf, 0), nrow = 2,
              dimnames = list(c("A", "B"), c("A", "B")))
  d <- stats::as.dist(m)
  expect_warning(
    res <- taxo_ordinate(d),
    "infinite values"
  )
  expect_null(res$points)
  expect_equal(res$dist, d)
  expect_s3_class(res, "taxodist_ord")
})

test_that("taxo_heatmap handles NA in distance matrix gracefully", {
  m <- matrix(c(0, NA, NA, 0), nrow = 2, dimnames = list(c("A", "B"), c("A", "B")))
  d <- stats::as.dist(m)
  pdf(file = NULL)
  expect_warning(res <- taxo_heatmap(d), "Distance matrix contains NA values")
  dev.off()
  expect_equal(res, d)
})

test_that("taxo_heatmap handles infinite distances gracefully", {
  m <- matrix(c(0, Inf, Inf, 0), nrow = 2,
              dimnames = list(c("A", "B"), c("A", "B")))
  d <- stats::as.dist(m)
  expect_warning(
    res <- taxo_heatmap(d),
    "infinite values"
  )
  expect_equal(res, d)
})

test_that("plot and summary methods safely ignore NULL components", {
  cl_null <- structure(list(hclust = NULL, dist = stats::dist(1:2)), class = "taxodist_cluster")
  expect_invisible(plot(cl_null))

  ord_null <- structure(list(points = NULL, dist = stats::dist(1:2), GOF = NULL, eig = NULL), class = "taxodist_ord")
  expect_invisible(plot(ord_null))
  expect_invisible(summary(ord_null))
})

