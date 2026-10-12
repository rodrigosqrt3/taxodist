# taxo_path

test_that("taxo_path returns NULL when taxon_a not found", {
  mockery::stub(taxo_path, "get_lineage", function(...) NULL)
  expect_null(taxo_path("Fakeosaurus", "Carnotaurus"))
})

test_that("taxo_path returns NULL when taxon_b not found", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Carnotaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Carnotaurus")
    else
      NULL
  })
  expect_null(taxo_path("Carnotaurus", "Fakeosaurus"))
})

test_that("taxo_path returns a taxodist_path data frame", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  })
  result <- taxo_path("Tyrannosaurus", "Triceratops")
  expect_s3_class(result, "taxodist_path")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("node", "depth", "direction"))
})

test_that("taxo_path has exactly one mrca row", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  })
  result <- taxo_path("Tyrannosaurus", "Triceratops")
  expect_equal(sum(result$direction == "mrca"), 1L)
})

test_that("taxo_path mrca node is correct", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  })
  result <- taxo_path("Tyrannosaurus", "Triceratops")
  expect_equal(result$node[result$direction == "mrca"], "Dinosauria")
})

test_that("taxo_path returns the complete path with lineage depths", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  })
  result <- taxo_path("Tyrannosaurus", "Triceratops")
  expect_equal(
    result$node,
    c("Tyrannosaurus", "Theropoda", "Dinosauria",
      "Ornithischia", "Triceratops")
  )
  expect_equal(result$depth, c(5L, 4L, 3L, 4L, 5L))
  expect_equal(result$direction, c("a", "a", "mrca", "b", "b"))
})

test_that("taxo_path returns only the MRCA for identical taxa", {
  lineage <- c("Biota", "Animalia", "Dinosauria", "Theropoda")
  mockery::stub(taxo_path, "get_lineage", function(...) lineage)
  result <- taxo_path("Theropoda", "Theropoda")
  expect_equal(result$node, "Theropoda")
  expect_equal(result$depth, 4L)
  expect_equal(result$direction, "mrca")
})

test_that("taxo_path handles taxon A as ancestor of taxon B", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Dinosauria")
      c("Biota", "Animalia", "Dinosauria")
    else
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
  })
  result <- taxo_path("Dinosauria", "Tyrannosaurus")
  expect_equal(result$node, c("Dinosauria", "Theropoda", "Tyrannosaurus"))
  expect_equal(result$depth, c(3L, 4L, 5L))
  expect_equal(result$direction, c("mrca", "b", "b"))
})

test_that("taxo_path handles taxon B as ancestor of taxon A", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria")
  })
  result <- taxo_path("Tyrannosaurus", "Dinosauria")
  expect_equal(result$node, c("Tyrannosaurus", "Theropoda", "Dinosauria"))
  expect_equal(result$depth, c(5L, 4L, 3L))
  expect_equal(result$direction, c("a", "a", "mrca"))
})

test_that("taxo_path returns NULL when no common ancestor exists", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Animalia") c("Biota", "Animalia")
    else c("Natura", "Mineralia")
  })
  expect_warning(
    result <- taxo_path("Animalia", "Mineralia"),
    "No common ancestor"
  )
  expect_null(result)
})

test_that("taxo_path direction column only contains valid values", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  })
  result <- taxo_path("Tyrannosaurus", "Triceratops")
  expect_true(all(result$direction %in% c("a", "mrca", "b")))
})

test_that("taxo_path preserves taxon attributes", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  })
  result <- taxo_path("Tyrannosaurus", "Triceratops")
  expect_equal(attr(result, "taxon_a"), "Tyrannosaurus")
  expect_equal(attr(result, "taxon_b"), "Triceratops")
})

