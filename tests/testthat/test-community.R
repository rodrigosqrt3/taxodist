# Taxonomic dissimilarity (Clarke, Somerfield & Chapman 2006) and
# standardized effect sizes of MPD and MNTD (Webb et al. 2002).

community_lineages <- list(
  A = c("Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "A"),
  B = c("Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "B"),
  C = c("Biota", "Animalia", "Chordata", "Mammalia", "Primates", "C"),
  D = c("Biota", "Animalia", "Arthropoda", "Insecta", "Diptera", "D"),
  E = c("Biota", "Animalia", "Arthropoda", "Insecta", "Coleoptera", "E"),
  F = c("Biota", "Plantae", "Tracheophyta", "Magnoliopsida", "Fagales", "F")
)
community_d <- function(metric = "path") {
  distance_matrix(taxo_from_lineages(community_lineages), metric = metric,
                  progress = FALSE)
}

test_that("Gamma+ and Theta+ reduce to Sorensen and Kulczynski on a flat classification", {
  skip_if_not_installed("vegan")
  set.seed(11)
  taxa <- paste0("sp", 1:12)
  flat <- stats::as.dist(1 - diag(12))
  attr(flat, "Labels") <- taxa
  comm <- matrix(rbinom(6 * 12, 1, 0.4), 6, dimnames = list(paste0("s", 1:6), taxa))
  comm[comm == 1] <- sample(1:9, sum(comm), replace = TRUE)  # abundances are ignored
  gamma <- taxo_dissimilarity(comm, flat, "gamma_plus")
  theta <- taxo_dissimilarity(comm, flat, "theta_plus")
  expect_equal(as.vector(gamma), as.vector(vegan::vegdist(comm > 0, "bray", binary = TRUE)))
  expect_equal(as.vector(theta), as.vector(vegan::vegdist(comm > 0, "kulczynski", binary = TRUE)))
  expect_identical(attr(gamma, "Labels"), rownames(comm))
  expect_identical(attr(gamma, "method"), "gamma_plus")
})

test_that("taxonomic dissimilarity follows its definition", {
  d <- community_d()
  comm <- rbind(
    s1 = c(A = 1, B = 0, C = 1, D = 0, E = 0, F = 0),
    s2 = c(A = 0, B = 1, C = 0, D = 0, E = 0, F = 0),
    s3 = c(A = 0, B = 0, C = 0, D = 1, E = 1, F = 0)
  )
  m <- as.matrix(d)
  # s1 vs s2: A -> B = 2, C -> B = 4; B -> A = 2.
  gamma <- as.matrix(taxo_dissimilarity(comm, d))
  theta <- as.matrix(taxo_dissimilarity(comm, d, "theta_plus"))
  expect_equal(gamma["s1", "s2"], (m["A", "B"] + m["C", "B"] + m["B", "A"]) / 3)
  expect_equal(theta["s1", "s2"], ((m["A", "B"] + m["C", "B"]) / 2 + m["B", "A"]) / 2)
  # Equal sample sizes give Gamma+ = Theta+.
  expect_equal(gamma["s1", "s3"], theta["s1", "s3"])
  expect_equal(attr(taxo_dissimilarity(comm, d), "metric"), "path")
  # Identical samples have zero dissimilarity.
  same <- rbind(a = c(A = 1, B = 1), b = c(A = 1, B = 1))
  expect_equal(as.vector(taxo_dissimilarity(same, d)), 0)
})

test_that("taxonomic dissimilarity handles empty samples and missing distances", {
  d <- community_d()
  comm <- rbind(s1 = c(A = 1, F = 0), s2 = c(A = 0, F = 0), s3 = c(A = 0, F = 1))
  values <- as.vector(taxo_dissimilarity(comm, d))
  # Pairs (s2, s1) and (s3, s2) involve the empty sample s2.
  expect_true(is.na(values[[1]]) && is.na(values[[3]]))
  expect_equal(values[[2]], as.matrix(d)["A", "F"])
  disconnected <- distance_matrix(taxo_from_lineages(list(
    X = c("Biota", "X"), Y = c("Other", "Y")
  )), progress = FALSE)
  comm2 <- rbind(a = c(X = 1, Y = 0), b = c(X = 0, Y = 1))
  expect_true(is.na(as.vector(taxo_dissimilarity(comm2, disconnected))))
  single <- taxo_dissimilarity(rbind(a = c(A = 1)), d)
  expect_length(single, 0L)
  expect_error(taxo_dissimilarity(comm, d, "jaccard"), "gamma_plus")
  expect_error(taxo_dissimilarity(character(0), d), "no taxa")
})

test_that("taxonomic dissimilarity computes distances when d is NULL", {
  calls <- 0L
  local_mocked_bindings(distance_matrix = function(taxa, ...) {
    calls <<- calls + 1L
    taxodist:::.distance_matrix_impl(
      taxo_from_lineages(community_lineages[taxa]), NULL,
      verbose = FALSE, progress = FALSE, metric = "path", levels = taxodist:::.rank_levels()
    )
  })
  out <- taxo_dissimilarity(rbind(a = c(A = 1, B = 0), b = c(A = 0, B = 1)), metric = "path")
  expect_identical(calls, 1L)
  expect_equal(as.vector(out), 2)
})

# Reference implementations written from picante's definitions.
reference_mpd <- function(x, m, weighted) {
  p <- x > 0
  w <- m[p, p]
  if (weighted) stats::weighted.mean(w, outer(x[p], x[p])) else mean(w[lower.tri(w)])
}
reference_mntd <- function(x, m, weighted) {
  p <- x > 0
  w <- m[p, p]
  diag(w) <- NA
  nearest <- apply(w, 2, min, na.rm = TRUE)
  if (weighted) stats::weighted.mean(nearest, x[p]) else mean(nearest)
}

test_that("MPD and MNTD follow picante's definitions", {
  d <- community_d()
  m <- as.matrix(d)
  comm <- rbind(
    s1 = c(A = 3, B = 1, C = 0, D = 2, E = 0, F = 1),
    s2 = c(A = 0, B = 5, C = 2, D = 0, E = 1, F = 0)
  )
  stat <- taxodist:::.community_statistic
  for (k in 1:2) {
    x <- comm[k, ]
    for (weighted in c(FALSE, TRUE)) {
      expect_equal(stat(x, m[names(x), names(x)], "mpd", weighted), reference_mpd(x, m[names(x), names(x)], weighted))
      expect_equal(stat(x, m[names(x), names(x)], "mntd", weighted), reference_mntd(x, m[names(x), names(x)], weighted))
    }
  }
  expect_true(is.na(stat(c(A = 1, B = 0), m[c("A", "B"), c("A", "B")], "mpd", FALSE)))
  inf <- m[c("A", "B"), c("A", "B")]
  inf[1, 2] <- inf[2, 1] <- Inf
  expect_true(is.na(stat(c(A = 1, B = 1), inf, "mpd", FALSE)))
})

test_that("the pool null model has the pool mean as expected MPD", {
  d <- community_d()
  set.seed(3)
  res <- taxo_ses(c("A", "B", "C"), d, n_iter = 4000)
  m <- as.matrix(d)
  expect_equal(res$null_mean, mean(m[lower.tri(m)]), tolerance = 0.02)
  expect_equal(res$observed, reference_mpd(c(A = 1, B = 1, C = 1), m[c("A", "B", "C"), c("A", "B", "C")], FALSE))
  expect_equal(res$NRI, -(res$observed - res$null_mean) / res$null_sd)
  # Three closely related mammals are clustered relative to the pool.
  expect_gt(res$NRI, 0)
  expect_s3_class(res, "taxodist_ses")
  expect_identical(attr(res, "pool_size"), 6L)
  expect_output(print(res), "sample")
  expect_message(print(res), "Net relatedness index")
})

test_that("taxo_ses supports MNTD, abundances, and custom pools", {
  d <- community_d()
  comm <- rbind(
    s1 = c(A = 3, B = 1, C = 0, D = 0),
    s2 = c(A = 0, B = 0, C = 2, D = 5),
    s3 = c(A = 1, B = 0, C = 0, D = 0)
  )
  set.seed(5)
  res <- taxo_ses(comm, d, index = "mntd", abundance_weighted = TRUE,
                  pool = c("A", "B", "C", "D", "E"), n_iter = 99)
  expect_named(res, c("sample", "S", "observed", "null_mean", "null_sd", "z",
                      "NTI", "p_lower", "p_upper", "n_iter"))
  expect_identical(res$S, c(2L, 2L, 1L))
  expect_true(is.na(res$observed[[3]]))
  expect_true(is.na(res$p_lower[[3]]))
  expect_true(all(res$p_lower[1:2] > 0 & res$p_lower[1:2] <= 1))
  expect_identical(attr(res, "pool_size"), 5L)
  expect_message(print(res), "abundance-weighted")
  expect_message(print(res), "Nearest taxon index")
})

test_that("the independent swap keeps richness and frequencies", {
  set.seed(9)
  comm <- matrix(rbinom(8 * 10, 1, 0.4), 8, dimnames = list(NULL, paste0("t", 1:10)))
  comm[comm == 1] <- sample(1:5, sum(comm), replace = TRUE)
  random <- taxodist:::.independent_swap(comm, 200)
  expect_identical(rowSums(random > 0), rowSums(comm > 0))
  expect_identical(colSums(random > 0), colSums(comm > 0))
  expect_false(identical(random > 0, comm > 0))
  # Abundances move with presences.
  expect_equal(sort(random[random > 0]), sort(comm[comm > 0]))
  # A matrix without checkerboards cannot change.
  full <- matrix(1, 2, 2, dimnames = list(NULL, c("a", "b")))
  expect_identical(taxodist:::.independent_swap(full, 5), full)
})

test_that("taxo_ses runs the independent swap null model", {
  d <- community_d()
  comm <- rbind(
    s1 = c(A = 1, B = 1, C = 0, D = 0, E = 0, F = 1),
    s2 = c(A = 0, B = 1, C = 1, D = 1, E = 0, F = 0),
    s3 = c(A = 1, B = 0, C = 0, D = 1, E = 1, F = 0)
  )
  set.seed(2)
  res <- taxo_ses(comm, d, null_model = "independent_swap", n_iter = 49, swaps = 50)
  expect_identical(nrow(res), 3L)
  expect_true(all(is.finite(res$null_mean)))
  expect_identical(attr(res, "null_model"), "independent_swap")
  expect_error(taxo_ses(c("A", "B"), d, null_model = "independent_swap"), "at least two samples")
})

test_that("taxo_ses validates its arguments", {
  d <- community_d()
  expect_error(taxo_ses(c("A", "B"), d, index = "pd"), "mpd")
  expect_error(taxo_ses(c("A", "B"), d, null_model = "frequency"), "pool")
  expect_error(taxo_ses(c("A", "B"), d, n_iter = 0), "positive integer")
  expect_error(taxo_ses(c("A", "B"), d, swaps = 1.5), "positive integer")
  expect_error(taxo_ses(c("A", "B"), d, abundance_weighted = NA), "TRUE or FALSE")
  expect_error(taxo_ses(c("A", "B"), d, pool = "A"), "at least two")
  expect_error(taxo_ses(c("A", "B"), d, pool = c("C", "D")), "belong to the pool")
  zero_sd <- taxo_ses(c("A", "B"), d, pool = c("A", "B"), n_iter = 5)
  expect_true(is.na(zero_sd$z))
  # A user-built dist object has no metric.
  plain <- stats::as.dist(as.matrix(d))
  attr(plain, "method") <- NULL
  res <- taxo_ses(c("A", "B", "C"), plain, n_iter = 5)
  expect_true(is.na(attr(res, "metric")))
  expect_false(any(grepl("metric", cli::cli_fmt(print(res)))))
})

test_that("taxo_ses computes distances when d is NULL", {
  local_mocked_bindings(distance_matrix = function(taxa, ...) {
    taxodist:::.distance_matrix_impl(
      taxo_from_lineages(community_lineages[taxa]), NULL,
      verbose = FALSE, progress = FALSE, metric = "reciprocal",
      levels = taxodist:::.rank_levels()
    )
  })
  set.seed(1)
  res <- taxo_ses(rbind(s = c(A = 1, B = 1, D = 1)), n_iter = 9)
  expect_identical(attr(res, "pool_size"), 3L)
  expect_identical(attr(res, "metric"), "reciprocal")
})

test_that("community indices reproduce the shared fixture", {
  fx <- jsonlite::read_json(test_path("fixtures", "expected_community.json"))
  lineages <- lapply(fx$lineages, unlist)
  res <- taxo_from_lineages(lineages)
  comm <- do.call(rbind, lapply(fx$comm, function(row) unlist(row)))
  rownames(comm) <- unlist(fx$sample_names)
  as_num <- function(x) vapply(x, function(v) if (is.null(v)) NA_real_ else v, 0)
  for (metric in c("reciprocal", "path")) {
    d <- distance_matrix(res, metric = metric, progress = FALSE)
    full <- as.matrix(d)
    exp <- fx[[metric]]
    expect_equal(as.vector(taxo_dissimilarity(comm, d, "gamma_plus")), as_num(exp$gamma_plus))
    expect_equal(as.vector(taxo_dissimilarity(comm, d, "theta_plus")), as_num(exp$theta_plus))
    for (index in c("mpd", "mntd")) for (weighted in c(FALSE, TRUE)) {
      key <- paste0(index, if (weighted) "_weighted" else "")
      got <- vapply(seq_len(nrow(comm)), function(k) {
        taxodist:::.community_statistic(comm[k, ], full, index, weighted)
      }, 0)
      expect_equal(got, as_num(exp[[key]]), info = paste(metric, key))
    }
  }
})
