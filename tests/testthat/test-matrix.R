# ── Vectorised distances, cross distances, and related helpers ──────────────

random_lineages <- function(n, depth = 8, branching = 3, seed = 1) {
  set.seed(seed)
  lineages <- lapply(seq_len(n), function(i) {
    len <- sample(2:depth, 1)
    path <- sample.int(branching, len, replace = TRUE)
    c("Biota", vapply(seq_len(len), function(d) {
      paste0("N", paste(path[seq_len(d)], collapse = ""))
    }, character(1)))
  })
  names(lineages) <- paste0("t", seq_len(n))
  lineages
}

test_that("vectorised distance_matrix matches pairwise .compute_distance", {
  lineages <- random_lineages(40)
  lineages$dup <- lineages$t1
  lineages$ancestor <- lineages$t2[1:2]
  lineages$foreign <- c("Natura", "Mineralia")
  matrix <- as.matrix(distance_matrix(taxo_from_lineages(lineages), progress = FALSE))

  for (i in seq_along(lineages)) {
    for (j in seq_along(lineages)) {
      expect_equal(
        matrix[i, j],
        taxodist:::.compute_distance(lineages[[i]], lineages[[j]])$distance
      )
    }
  }
})

test_that("prefix node IDs handle empty input and lineages without nodes", {
  expect_equal(dim(taxodist:::.prefix_node_ids(list())), c(0L, 0L))
  expect_equal(
    taxodist:::.pairwise_distances(list(character(0), character(0))),
    Inf
  )
  expect_equal(taxodist:::.pairwise_distances(list(NULL, NULL)), NA_real_)
})

test_that("distance_matrix rejects unsupported inputs and tolerates NA names", {
  expect_error(distance_matrix(list("a", "b")), "character vector")
  mockery::stub(distance_matrix, "get_lineage", function(taxon, ...) {
    c("Biota", taxon)
  })
  result <- as.matrix(distance_matrix(c("Alpha", NA, "Beta"), progress = FALSE))
  expect_true(is.na(result[1, 2]))
  expect_equal(result[1, 3], 1)
})

test_that("distance_matrix accepts factors and reuses its name cache", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  calls <- 0L
  mockery::stub(distance_matrix, "get_lineage", function(taxon, ...) {
    calls <<- calls + 1L
    c("Biota", taxon)
  })
  first <- distance_matrix(factor(c("Alpha", "Beta")), progress = FALSE)
  second <- distance_matrix(c(" Alpha", "Beta "), progress = FALSE)
  expect_equal(calls, 2L)
  expect_equal(attr(first, "Labels"), c("Alpha", "Beta"))
  expect_equal(c(first), c(second))
})

test_that("cross_distances returns the rectangular block of distance_matrix", {
  lineages <- random_lineages(12, seed = 7)
  x <- taxo_from_lineages(lineages[1:5])
  y <- taxo_from_lineages(lineages[6:12])
  full <- as.matrix(distance_matrix(taxo_from_lineages(lineages), progress = FALSE))

  cross <- cross_distances(x, y, progress = FALSE)
  expect_equal(dim(cross), c(5L, 7L))
  expect_equal(rownames(cross), names(lineages)[1:5])
  expect_equal(colnames(cross), names(lineages)[6:12])
  expect_equal(cross, full[1:5, 6:12])
})

test_that("cross_distances marks unresolved taxa as NA and handles empty sets", {
  x <- data.frame(
    input = c("Alpha", "Missing"),
    status = c("resolved", "unresolved"),
    stringsAsFactors = FALSE
  )
  x$lineage <- I(list(c("Biota", "Alpha"), NULL))
  class(x) <- c("taxodist_resolution", "data.frame")
  y <- x
  y$input <- c("Beta", "Gone")
  y$lineage <- I(list(c("Biota", "Beta"), NULL))

  result <- cross_distances(x, y, progress = FALSE)
  expect_equal(result[1, 1], 1)
  expect_true(all(is.na(result[2, ])))
  expect_true(all(is.na(result[, 2])))

  empty <- cross_distances(character(0), x, progress = FALSE)
  expect_equal(dim(empty), c(0L, 2L))
})

test_that("cross_distances retrieves character inputs", {
  mockery::stub(cross_distances, "get_lineage", function(taxon, ...) {
    c("Biota", "Animalia", taxon)
  })
  result <- cross_distances("Alpha", c("Beta", "Alpha"), progress = TRUE)
  expect_equal(unname(result[1, ]), c(0.5, 0))
})

test_that("bundles can be used wherever a resolution is accepted", {
  bundle <- taxo_bundle(taxo_from_lineages(list(
    Alpha = c("Biota", "Alpha"),
    Beta = c("Biota", "Beta")
  )), progress = FALSE)
  expect_equal(cross_distances(bundle, bundle, progress = FALSE)[1, 2], 1)
})

test_that("closest_relative returns consecutive row names", {
  mockery::stub(closest_relative, "get_lineage", function(taxon, ...) {
    list(
      Query = c("Biota", "Animalia", "Chordata", "Query"),
      Far = c("Biota", "Plantae", "Far"),
      Near = c("Biota", "Animalia", "Chordata", "Near")
    )[[taxon]]
  })
  result <- closest_relative("Query", c("Far", "Near"))
  expect_equal(result$taxon, c("Near", "Far"))
  expect_equal(rownames(result), c("1", "2"))
})

test_that("focal_distances reports the lineage node as MRCA for identical taxa", {
  mockery::stub(focal_distances, "get_lineage", function(...) {
    c("Biota", "Dinosauria", "Tyrannosaurus")
  })
  result <- focal_distances("tyrannosaurus", "Tyrannosaurus", progress = FALSE)
  expect_equal(result$distance, 0)
  expect_equal(result$mrca, "Tyrannosaurus")
  expect_equal(result$mrca_depth, 3L)
})

test_that("taxo_heatmap passes graphical arguments only to heatmap", {
  mockery::stub(taxo_heatmap, "distance_matrix", function(taxa, verbose, progress,
                                                          metric = "reciprocal") {
    expect_false(progress)
    distance_matrix(taxo_from_lineages(list(
      A = c("Biota", "X", "A"),
      B = c("Biota", "X", "B"),
      C = c("Biota", "Y", "C")
    )), progress = FALSE)
  })
  pdf(file = NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_warning(
    result <- taxo_heatmap(
      c("A", "B", "C"),
      col = grDevices::heat.colors(5),
      main = "Test",
      progress = FALSE
    )
  )
  expect_s3_class(result, "dist")
})

test_that("sorted-prefix MRCA index agrees with pairwise comparison", {
  set.seed(20261124)
  random_lineage <- function() {
    depth <- sample(1:12, 1)
    paste0("n", seq_len(depth), "_", sample(0:1, depth, TRUE))
  }
  lineages <- c(
    replicate(40, random_lineage(), simplify = FALSE),
    list(
      NULL,                                   # unresolved taxon
      character(0),                           # empty lineage
      c("n1_0", "n2_0"),                      # ancestor of others
      c("n1_0", "n2_0"),                      # duplicate
      c("Other", "root"),                     # disconnected hierarchy
      c("n1_0", "n2_0", "n3_1", "n4_0", "n5_1")
    )
  )
  n <- length(lineages)
  for (metric in c("reciprocal", "path", "wu_palmer")) {
    fast <- taxodist:::.pairwise_distances(lineages, metric = metric)
    slow <- numeric(0)
    for (j in seq_len(n - 1L)) for (i in (j + 1L):n) {
      a <- lineages[[i]]
      b <- lineages[[j]]
      slow <- c(slow, if (is.null(a) || is.null(b)) NA_real_ else
        taxodist:::.compute_distance(a, b, metric = metric)$distance)
    }
    expect_equal(fast, slow, info = metric)
  }

  index <- taxodist:::.mrca_index(taxodist:::.prefix_node_ids(lineages))
  for (j in c(1L, 20L, n - 1L, n)) {
    depths <- taxodist:::.mrca_depths_from(index, j)
    expected <- vapply(lineages, function(l) {
      if (is.null(l) || is.null(lineages[[j]])) 0L else
        taxodist:::.common_prefix_depth(l, lineages[[j]])
    }, integer(1))
    expect_identical(depths, expected, info = j)
  }
})

test_that("cross_distances uses the same index as distance_matrix", {
  set.seed(7)
  lineages <- replicate(25, paste0("c", seq_len(sample(2:8, 1)), sample(0:1, 8, TRUE)[1:8])[1:sample(2:8, 1)],
                        simplify = FALSE)
  names(lineages) <- paste0("t", seq_along(lineages))
  res <- taxo_from_lineages(lineages)
  full <- as.matrix(distance_matrix(res, progress = FALSE, metric = "wu_palmer"))
  x <- taxo_from_lineages(lineages[1:10])
  y <- taxo_from_lineages(lineages[11:25])
  cross <- cross_distances(x, y, progress = FALSE, metric = "wu_palmer")
  expect_equal(unname(cross), unname(full[1:10, 11:25]))
})

test_that("MRCA index handles degenerate inputs", {
  empty <- taxodist:::.mrca_index(matrix(integer(0), nrow = 0, ncol = 0))
  expect_identical(empty$n, 0L)
  single <- taxodist:::.mrca_index(taxodist:::.prefix_node_ids(list(c("a", "b"))))
  expect_identical(taxodist:::.mrca_depths_from(single, 1L), 2L)
  flat <- taxodist:::.mrca_index(taxodist:::.prefix_node_ids(list(character(0), character(0))))
  expect_identical(taxodist:::.mrca_depths_from(flat, 1L), c(0L, 0L))
})
