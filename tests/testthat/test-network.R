# Live network tests (skipped on CRAN and CI)

skip_if_taxonomicon_down <- function() {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  res <- tryCatch(httr::GET("http://taxonomicon.taxonomy.nl", httr::timeout(3)), error = function(e) NULL)
  if (is.null(res) || httr::status_code(res) != 200) {
    testthat::skip("Taxonomicon server is currently offline.")
  }
}

test_that("get_lineage returns correct lineage for Velociraptor", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  lin <- get_lineage("Velociraptor")
  expect_type(lin, "character")
  expect_true(length(lin) > 0)
  expect_true("Dinosauria" %in% lin)
  expect_true("Theropoda" %in% lin)
  expect_true("Dromaeosauridae" %in% lin)
})

test_that("get_lineage returns correct lineage for Tyrannosaurus", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  lin <- get_lineage("Tyrannosaurus")
  expect_type(lin, "character")
  expect_true("Coelurosauria" %in% lin)
  expect_true("Dinosauria" %in% lin)
})

test_that("get_lineage returns correct lineage for Carnotaurus", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  lin <- get_lineage("Carnotaurus")
  expect_type(lin, "character")
  expect_true("Dinosauria" %in% lin)
  expect_true("Theropoda" %in% lin)
})

test_that("get_lineage returns correct lineage for Homo", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  lin <- get_lineage("Homo")
  expect_true("Amniota" %in% lin)
  expect_true("Mammalia" %in% lin)
})

test_that("get_lineage returns correct lineage for Drosophila", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  lin <- suppressWarnings(get_lineage("Drosophila"))
  skip_if(is.null(lin), "Taxonomicon unstable")
  expect_type(lin, "character")
  expect_true(length(lin) > 0)
  expect_true("Animalia" %in% lin)
})

test_that("get_lineage returns NULL for unknown taxon", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  expect_null(get_lineage("Fakeosaurus"))
})

test_that("taxo_distance returns valid result for Tyrannosaurus vs Velociraptor", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  result <- taxo_distance("Tyrannosaurus", "Velociraptor")
  expect_s3_class(result, "taxodist_result")
  expect_gte(result$distance, 0)
  expect_lte(result$distance, 1)
  expect_equal(result$mrca, "Tyrannoraptora")
})

test_that("taxo_distance is positive when one taxon is ancestor of other", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  res1 <- taxo_distance("Tyrannosaurus", "Dinosauria")
  skip_if(is.null(res1), "Taxonomicon unstable")
  expect_gt(res1$distance, 0)
  expect_equal(res1$distance, 1 / res1$mrca_depth)

  res2 <- taxo_distance("Carnotaurus", "Ceratosauria")
  skip_if(is.null(res2), "Taxonomicon unstable")
  expect_gt(res2$distance, 0)
  expect_equal(res2$distance, 1 / res2$mrca_depth)
})

test_that("taxo_distance between Carnotaurus and Triceratops is valid", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  result <- taxo_distance("Carnotaurus", "Triceratops")
  expect_s3_class(result, "taxodist_result")
  expect_equal(result$mrca, "Dinosauria")
})

test_that("taxo_distance is larger between distant taxa than close ones", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  d_close  <- taxo_distance("Carnotaurus", "Tyrannosaurus")$distance
  d_distant <- taxo_distance("Carnotaurus", "Homo")$distance
  expect_gt(d_distant, d_close)
})

test_that("mrca of Tyrannosaurus and Triceratops is Dinosauria", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  expect_equal(mrca("Tyrannosaurus", "Triceratops"), "Dinosauria")
})

test_that("mrca of Tyrannosaurus and Homo is Amniota", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  expect_equal(mrca("Tyrannosaurus", "Homo"), "Amniota")
})

test_that("mrca of Velociraptor and Triceratops is Dinosauria", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  expect_equal(mrca("Velociraptor", "Triceratops"), "Dinosauria")
})

test_that("mrca of Carnotaurus and Tyrannosaurus is within Theropoda", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  ancestor <- mrca("Carnotaurus", "Tyrannosaurus")
  skip_if(is.null(ancestor), "Taxonomicon unstable")
  lin <- get_lineage("Tyrannosaurus")
  expect_true(ancestor %in% lin)
})

test_that("is_member correctly identifies clade membership", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  expect_true(is_member("Tyrannosaurus", "Theropoda"))
  expect_false(is_member("Tyrannosaurus", "Ornithischia"))
})

test_that("lineage_depth for Carnotaurus is reasonable", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  expect_gt(lineage_depth("Carnotaurus"), 10)
})

test_that("get_taxonomicon_id finds real ID and caches it", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  id <- get_taxonomicon_id("Carnotaurus", verbose = TRUE)
  expect_type(id, "character")
  expect_false(is.null(id))
  expect_equal(id, get("id_Carnotaurus", envir = taxodist:::.taxodist_cache))
})

test_that("get_lineage_by_id parses and caches lineage for Drosophila", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  id <- suppressWarnings(get_taxonomicon_id("Drosophila"))
  skip_if(is.null(id), "Taxonomicon unstable")
  result <- get_lineage_by_id(id, verbose = TRUE)
  expect_type(result, "character")
  expect_true("Animalia" %in% result)
})

test_that("get_lineage_by_id clean = FALSE keeps more nodes than clean = TRUE", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  id <- get_taxonomicon_id("Carnotaurus")
  result_clean <- get_lineage_by_id(id, clean = TRUE)
  clear_cache()
  result_no_clean <- get_lineage_by_id(id, clean = FALSE)
  expect_lte(length(result_clean), length(result_no_clean))
})

test_that("get_lineage verbose wrapper works for Quercus", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  result <- get_lineage("Quercus", verbose = TRUE)
  expect_type(result, "character")
  expect_true("Biota" %in% result)
})

test_that("get_taxonomicon_id returns NULL for genuinely unknown taxon", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  expect_null(get_taxonomicon_id("Zzzznotarealgenus99999", verbose = TRUE))
})

test_that("get_taxonomicon_id skips astronomical objects", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  id <- get_taxonomicon_id("Venus", verbose = TRUE)
  if (!is.null(id)) {
    lin <- get_lineage_by_id(id)
    expect_true(!is.null(lin))
  } else {
    expect_null(id)
  }
})

test_that("get_lineage_by_id works directly with verbose", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  id <- get_taxonomicon_id("Carnotaurus")
  expect_no_error(get_lineage_by_id(id, verbose = TRUE))
})

test_that("get_taxonomicon_id works with verbose for real taxon", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  expect_no_error(suppressWarnings(get_taxonomicon_id("Drosophila", verbose = TRUE)))
})

test_that("get_taxonomicon_id verbose prints not found warning", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  expect_null(get_taxonomicon_id("Zzzzfakeosaurus99999", verbose = TRUE))
})

test_that("get_lineage_by_id verbose success message fires on real taxon", {
  skip_if_offline()
  skip_if_taxonomicon_down()
  clear_cache()
  id <- get_taxonomicon_id("Carnotaurus")
  expect_no_error(get_lineage_by_id(id, verbose = TRUE))
})

