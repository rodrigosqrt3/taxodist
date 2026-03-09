test_that("taxodist package loads", {
  expect_true(TRUE)
})

test_that(".compute_distance works correctly", {
  lin_a <- c("Biota", "Animalia", "Chordata", "Dinosauria",
             "Theropoda", "Tyrannosauridae", "Tyrannosaurus")
  lin_b <- c("Biota", "Animalia", "Chordata", "Dinosauria",
             "Theropoda", "Dromaeosauridae", "Velociraptor")

  result <- taxodist:::.compute_distance(lin_a, lin_b, "Tyrannosaurus",
                                          "Velociraptor", method = "raw")

  expect_equal(result$mrca, "Theropoda")
  expect_equal(result$mrca_depth, 5L)
  expect_equal(result$distance, 4L)  # 7 + 7 - 2*5
  expect_equal(result$depth_a, 7L)
  expect_equal(result$depth_b, 7L)
})

test_that(".compute_distance jaccard is between 0 and 1", {
  lin_a <- c("Biota", "Animalia", "Chordata", "Dinosauria",
             "Theropoda", "Tyrannosauridae", "Tyrannosaurus")
  lin_b <- c("Biota", "Animalia", "Chordata", "Dinosauria",
             "Theropoda", "Dromaeosauridae", "Velociraptor")

  result <- taxodist:::.compute_distance(lin_a, lin_b, method = "jaccard")
  expect_gte(result$distance, 0)
  expect_lte(result$distance, 1)
})

test_that(".compute_distance is symmetric", {
  lin_a <- c("Biota", "Animalia", "Chordata", "Dinosauria", "Theropoda")
  lin_b <- c("Biota", "Animalia", "Chordata", "Dinosauria", "Ornithischia")

  r1 <- taxodist:::.compute_distance(lin_a, lin_b, method = "raw")
  r2 <- taxodist:::.compute_distance(lin_b, lin_a, method = "raw")

  expect_equal(r1$distance, r2$distance)
})

test_that(".compute_distance satisfies triangle inequality (raw)", {
  lin_a <- c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosauridae")
  lin_b <- c("Biota", "Animalia", "Dinosauria", "Theropoda", "Dromaeosauridae")
  lin_c <- c("Biota", "Animalia", "Dinosauria", "Ornithischia")

  dAB <- taxodist:::.compute_distance(lin_a, lin_b, method = "raw")$distance
  dBC <- taxodist:::.compute_distance(lin_b, lin_c, method = "raw")$distance
  dAC <- taxodist:::.compute_distance(lin_a, lin_c, method = "raw")$distance

  expect_lte(dAC, dAB + dBC)
})

test_that(".compute_distance returns 0 for identical lineages", {
  lin <- c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus")
  result <- taxodist:::.compute_distance(lin, lin, method = "raw")
  expect_equal(result$distance, 0L)
  expect_equal(result$mrca, "Tyrannosaurus")
})

test_that(".compute_distance handles no common ancestor", {
  lin_a <- c("Biota", "Animalia")
  lin_b <- c("Fungi", "Ascomycota")
  result <- taxodist:::.compute_distance(lin_a, lin_b, method = "raw")
  expect_equal(result$mrca_depth, 0L)
  expect_true(is.na(result$mrca))
})

test_that("filter_clade filters correctly with mock lineages", {
  # mock get_lineage to avoid network calls in tests
  mockery::stub(filter_clade, "is_member", function(taxon, clade, ...) {
    memberships <- list(
      Tyrannosaurus = c("Dinosauria", "Theropoda"),
      Triceratops   = c("Dinosauria", "Ornithischia"),
      Homo          = c("Mammalia", "Amniota")
    )
    clade %in% memberships[[taxon]]
  })

  result <- filter_clade(
    c("Tyrannosaurus", "Triceratops", "Homo"), "Dinosauria"
  )
  expect_equal(result, c("Tyrannosaurus", "Triceratops"))
})
