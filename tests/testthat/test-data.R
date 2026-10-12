# Package data

test_that("taxodist package loads", {
  expect_true(TRUE)
})

test_that("taxobase has the documented structure", {
  data("taxobase", package = "taxodist", envir = environment())
  expect_type(taxobase, "list")
  expect_named(taxobase, c(
    "taxa", "found_taxa", "coverage", "matrix", "pairwise",
    "lineage_homo", "lineage_tyrannosaurus", "closest", "filter",
    "search", "statistical_taxa", "statistical_matrix", "metadata"
  ))
  expect_s3_class(taxobase$matrix, "dist")
  expect_s3_class(taxobase$statistical_matrix, "dist")
  expect_equal(attr(taxobase$matrix, "Labels"), taxobase$found_taxa)
  expect_equal(
    attr(taxobase$statistical_matrix, "Labels"),
    taxobase$statistical_taxa
  )
  expect_type(taxobase$metadata$package_version, "character")
  expect_length(taxobase$metadata$package_version, 1L)
  expect_true(nzchar(taxobase$metadata$package_version))
  expect_true(
    base::package_version(taxobase$metadata$package_version) <=
      utils::packageVersion("taxodist")
  )
})

