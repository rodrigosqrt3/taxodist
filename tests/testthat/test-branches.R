# Alternative branches of inline conditionals. covr >= 4.0 counts every
# branch of `if (...) a else b` written on one line, so each one is
# exercised here explicitly.

test_that("metric specs that are not lists are rejected", {
  expect_error(taxodist:::.metric_from_spec("inverse_mrca_depth"), "unknown distance metric")
  expect_error(taxodist:::.metric_from_spec(NULL), "unknown distance metric")
})

test_that("bundle matrices encode both infinities", {
  m <- structure(
    c(Inf, -Inf, NA),
    Size = 3L, Labels = c("a", "b", "c"), Diag = FALSE, Upper = FALSE,
    class = "dist"
  )
  rows <- taxodist:::.bundle_matrix_rows(m)
  expect_identical(rows[[2]][[1]], "Infinity")
  expect_identical(rows[[3]][[1]], "-Infinity")
  expect_null(rows[[3]][[2]])
})

test_that("distinctness test prints without a metric", {
  pool <- taxo_from_lineages(
    list(
      A = c("Biota", "Animalia", "Chordata", "A"),
      B = c("Biota", "Animalia", "Chordata", "B"),
      C = c("Biota", "Animalia", "Arthropoda", "C"),
      D = c("Biota", "Plantae", "D")
    )
  )
  d <- distance_matrix(pool, progress = FALSE)
  set.seed(1)
  result <- taxo_distinctness_test(c("A", "C"), d, n_iter = 19)
  result$metric <- NA_character_
  out <- cli::cli_fmt(print(result))
  expect_false(any(grepl("metric:", out, fixed = TRUE)))
  result$metric <- "reciprocal"
  out <- cli::cli_fmt(print(result))
  expect_true(any(grepl("metric: reciprocal", out, fixed = TRUE)))
})

test_that("rank text without the name falls back to the form of the name", {
  expect_identical(taxodist:::.rank_before_name("Family Felidae", "Felidae"), "Family")
  expect_identical(taxodist:::.rank_before_name("nothing here", "Homo sapiens"), "Species")
})

test_that("raw ranks are cached under their own key", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  fake <- function(taxon_id, clean = TRUE, verbose = FALSE) {
    lineage <- c("Natura", "Biota", "Alpha")
    assign(paste0(if (clean) "ranks_" else "rawranks_", taxon_id),
           c(NA, NA, "Genus"), envir = taxodist:::.taxodist_cache)
    lineage
  }
  ranks_by_id <- taxodist:::.ranks_by_id
  mockery::stub(ranks_by_id, "get_lineage_by_id", fake)
  expect_identical(ranks_by_id("123", clean = FALSE), c(NA, NA, "Genus"))
  expect_true(exists("rawranks_123", envir = taxodist:::.taxodist_cache))
})

test_that("results without a stored metric print as reciprocal", {
  res <- structure(
    list(taxon_a = "A", taxon_b = "B", distance = 0.5, mrca = "Biota",
         mrca_depth = 2L, depth_a = 3L, depth_b = 3L),
    class = "taxodist_result"
  )
  out <- cli::cli_fmt(print(res))
  expect_true(any(grepl("(reciprocal)", out, fixed = TRUE)))
})

test_that("empty cluster and ordination objects plot without error", {
  cl <- structure(list(hclust = NULL), class = "taxodist_cluster")
  expect_invisible(plot(cl))
  ord <- structure(list(points = NULL), class = "taxodist_ord")
  expect_invisible(plot(ord))
})
