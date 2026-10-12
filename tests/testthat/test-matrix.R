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
  mockery::stub(taxo_heatmap, "distance_matrix", function(taxa, verbose, progress) {
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
