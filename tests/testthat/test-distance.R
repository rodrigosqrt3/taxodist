library(testthat)
library(taxodist)
library(webmockr)

# ── Pure logic tests ──────────────────────────────────────────────────────────

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

test_that(".compute_distance works correctly", {
  lin_a <- c("Biota", "Animalia", "Chordata", "Dinosauria",
             "Theropoda", "Tyrannosauridae", "Tyrannosaurus")
  lin_b <- c("Biota", "Animalia", "Chordata", "Dinosauria",
             "Theropoda", "Dromaeosauridae", "Velociraptor")
  result <- taxodist:::.compute_distance(lin_a, lin_b, "Tyrannosaurus",
                                         "Velociraptor")
  expect_equal(result$mrca, "Theropoda")
  expect_equal(result$mrca_depth, 5L)
  expect_equal(result$depth_a, 7L)
  expect_equal(result$depth_b, 7L)
  expect_gte(result$distance, 0)
  expect_lte(result$distance, 1)
})

test_that(".compute_distance distance is between 0 and 1", {
  lin_a <- c("Biota", "Animalia", "Chordata", "Dinosauria",
             "Theropoda", "Tyrannosauridae", "Tyrannosaurus")
  lin_b <- c("Biota", "Animalia", "Chordata", "Dinosauria",
             "Theropoda", "Dromaeosauridae", "Velociraptor")
  result <- taxodist:::.compute_distance(lin_a, lin_b)
  expect_gte(result$distance, 0)
  expect_lte(result$distance, 1)
})

test_that(".compute_distance is symmetric", {
  lin_a <- c("Biota", "Animalia", "Chordata", "Dinosauria", "Theropoda")
  lin_b <- c("Biota", "Animalia", "Chordata", "Dinosauria", "Ornithischia")
  r1 <- taxodist:::.compute_distance(lin_a, lin_b)
  r2 <- taxodist:::.compute_distance(lin_b, lin_a)
  expect_equal(r1$distance, r2$distance)
})

test_that(".compute_distance ignores names repeated after divergence", {
  lin_a <- c("Biota", "Animalia", "Metazoa", "RepeatedName", "TaxonA")
  lin_b <- c("Biota", "Plantae", "Viridiplantae", "RepeatedName", "TaxonB")
  result <- taxodist:::.compute_distance(lin_a, lin_b)
  expect_equal(result$mrca, "Biota")
  expect_equal(result$mrca_depth, 1L)
  expect_equal(result$distance, 1)
})

test_that(".compute_distance uses only the continuous common prefix", {
  lin_a <- c("Biota", "Animalia", "Chordata", "SharedAgain", "TaxonA")
  lin_b <- c("Biota", "Animalia", "Arthropoda", "Insecta", "SharedAgain",
             "TaxonB")
  result <- taxodist:::.compute_distance(lin_a, lin_b)
  expect_equal(result$mrca, "Animalia")
  expect_equal(result$mrca_depth, 2L)
  expect_equal(result$distance, 0.5)
})

test_that(".compute_distance remains symmetric for unequal lineage lengths", {
  lin_a <- c("Biota", "Animalia", "Chordata", "SharedAgain", "TaxonA")
  lin_b <- c("Biota", "Animalia", "Arthropoda", "SharedAgain", "Group",
             "TaxonB")
  r1 <- taxodist:::.compute_distance(lin_a, lin_b)
  r2 <- taxodist:::.compute_distance(lin_b, lin_a)
  expect_equal(r1$distance, r2$distance)
  expect_equal(r1$mrca, r2$mrca)
  expect_equal(r1$mrca_depth, r2$mrca_depth)
})

test_that(".compute_distance satisfies triangle inequality", {
  lin_a <- c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosauridae")
  lin_b <- c("Biota", "Animalia", "Dinosauria", "Theropoda", "Dromaeosauridae")
  lin_c <- c("Biota", "Animalia", "Dinosauria", "Ornithischia")
  dAB <- taxodist:::.compute_distance(lin_a, lin_b)$distance
  dBC <- taxodist:::.compute_distance(lin_b, lin_c)$distance
  dAC <- taxodist:::.compute_distance(lin_a, lin_c)$distance
  expect_lte(dAC, dAB + dBC)
})

test_that(".compute_distance returns 0 for identical lineages", {
  lin <- c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus")
  result <- taxodist:::.compute_distance(lin, lin)
  expect_equal(result$distance, 0)
  expect_equal(result$mrca, "Tyrannosaurus")
})

test_that(".compute_distance handles no common ancestor", {
  lin_a <- c("Biota", "Animalia")
  lin_b <- c("Fungi", "Ascomycota")
  result <- taxodist:::.compute_distance(lin_a, lin_b)
  expect_equal(result$mrca_depth, 0L)
  expect_true(is.na(result$mrca))
})

test_that(".compute_distance returns Inf for no shared ancestor", {
  lin_a <- c("Biota", "Animalia")
  lin_b <- c("Fungi", "Ascomycota")
  result <- taxodist:::.compute_distance(lin_a, lin_b)
  expect_equal(result$distance, Inf)
})

test_that(".compute_distance result has correct S3 class", {
  lin <- c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus")
  result <- taxodist:::.compute_distance(lin, lin)
  expect_s3_class(result, "taxodist_result")
})

test_that(".compute_distance distance is between 0 and 1 for asymmetric lineages", {
  lin_a <- c("Biota", "Animalia", "Dinosauria", "Theropoda",
             "Abelisauridae", "Carnotaurus")
  lin_b <- c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  result <- taxodist:::.compute_distance(lin_a, lin_b)
  expect_gte(result$distance, 0)
  expect_lte(result$distance, 1)
})

test_that(".compute_distance is positive when one taxon is ancestor of other", {
  lin_a <- c("Biota", "Animalia", "Dinosauria")
  lin_b <- c("Biota", "Animalia", "Dinosauria", "Theropoda", "Carnotaurus")
  result <- taxodist:::.compute_distance(lin_a, lin_b)
  reverse <- taxodist:::.compute_distance(lin_b, lin_a)
  expect_equal(result$distance, 1 / 3)
  expect_equal(reverse$distance, result$distance)
  expect_equal(result$mrca, "Dinosauria")
})

test_that(".compute_distance satisfies ultrametric inequality with ancestor nodes", {
  ancestor <- c("Biota", "Animalia", "Dinosauria")
  theropod <- c("Biota", "Animalia", "Dinosauria", "Theropoda")
  ornithischian <- c("Biota", "Animalia", "Dinosauria", "Ornithischia")

  d_ab <- taxodist:::.compute_distance(ancestor, theropod)$distance
  d_ac <- taxodist:::.compute_distance(ancestor, ornithischian)$distance
  d_bc <- taxodist:::.compute_distance(theropod, ornithischian)$distance

  expect_gt(d_ab, 0)
  expect_lte(d_ab, max(d_ac, d_bc))
  expect_lte(d_ac, max(d_ab, d_bc))
  expect_lte(d_bc, max(d_ab, d_ac))
})

test_that("numeric distance path matches detailed distance results", {
  cases <- list(
    list(
      c("Biota", "Animalia", "Dinosauria"),
      c("Biota", "Animalia", "Dinosauria")
    ),
    list(
      c("Biota", "Animalia", "Dinosauria"),
      c("Biota", "Animalia", "Dinosauria", "Theropoda")
    ),
    list(
      c("Biota", "Animalia", "Chordata", "TaxonA"),
      c("Biota", "Animalia", "Arthropoda", "TaxonB")
    ),
    list(
      c("Biota", "Animalia", "Repeated", "TaxonA"),
      c("Biota", "Plantae", "Repeated", "TaxonB")
    ),
    list(
      c("Biota", "Animalia"),
      c("Natura", "Mineralia")
    )
  )

  for (case in cases) {
    expect_equal(
      taxodist:::.distance_value(case[[1L]], case[[2L]]),
      taxodist:::.compute_distance(case[[1L]], case[[2L]])$distance
    )
  }
})

test_that("common-prefix helper handles an empty lineage", {
  expect_equal(
    taxodist:::.common_prefix_depth(character(0), "Biota"),
    0L
  )
})

test_that("clear_cache returns invisible NULL", {
  expect_invisible(clear_cache())
})

test_that("save_cache creates a file with cache contents", {
  clear_cache()
  assign("id_Carnotaurus", "12345", envir = taxodist:::.taxodist_cache)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  expect_invisible(save_cache(tmp))
  expect_true(file.exists(tmp))
})

test_that("load_cache restores entries into the cache", {
  clear_cache()
  assign("id_Carnotaurus", "12345", envir = taxodist:::.taxodist_cache)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  save_cache(tmp)
  clear_cache()
  expect_false(exists("id_Carnotaurus", envir = taxodist:::.taxodist_cache))
  load_cache(tmp)
  expect_true(exists("id_Carnotaurus", envir = taxodist:::.taxodist_cache))
  expect_equal(get("id_Carnotaurus", envir = taxodist:::.taxodist_cache), "12345")
})

test_that("save_cache / load_cache round-trip preserves all entries", {
  clear_cache()
  assign("id_Tyrannosaurus", "50841", envir = taxodist:::.taxodist_cache)
  assign("lin_50841", c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  save_cache(tmp)
  clear_cache()
  load_cache(tmp)
  expect_equal(get("id_Tyrannosaurus",  envir = taxodist:::.taxodist_cache), "50841")
  expect_equal(get("lin_50841", envir = taxodist:::.taxodist_cache),
               c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus"))
})

test_that("load_cache errors on missing file", {
  expect_error(load_cache("nonexistent_file.rds"))
})

test_that("load_cache rejects a non-list RDS without modifying cache", {
  clear_cache()
  assign("id_existing", "123", envir = taxodist:::.taxodist_cache)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(42, tmp)
  expect_error(load_cache(tmp), "expected a named list")
  expect_equal(
    get("id_existing", envir = taxodist:::.taxodist_cache),
    "123"
  )
})

test_that("load_cache rejects unnamed cache entries", {
  clear_cache()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(list("12345"), tmp)
  expect_error(load_cache(tmp), "unique, non-empty names")
})

test_that("load_cache rejects invalid ID and lineage values", {
  clear_cache()
  tmp_id <- tempfile(fileext = ".rds")
  tmp_lin <- tempfile(fileext = ".rds")
  tmp_matrix_lin <- tempfile(fileext = ".rds")
  on.exit(unlink(c(tmp_id, tmp_lin, tmp_matrix_lin)))

  saveRDS(list(id_bad = 12345), tmp_id)
  expect_error(load_cache(tmp_id), "ID entries")

  saveRDS(list(lin_bad = 12345), tmp_lin)
  expect_error(load_cache(tmp_lin), "lineage entries")

  saveRDS(list(matrix_lineage_bad = 12345), tmp_matrix_lin)
  expect_error(load_cache(tmp_matrix_lin), "matrix lineage entries")
})

test_that("save_cache returns invisible NULL", {
  clear_cache()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  expect_invisible(save_cache(tmp))
})

test_that("load_cache returns invisible NULL", {
  clear_cache()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  save_cache(tmp)
  clear_cache()
  expect_invisible(load_cache(tmp))
  expect_equal(ls(taxodist:::.taxodist_cache), character(0))
})

test_that("filter_clade filters correctly with mocked lineages", {
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

# ── taxo_path ─────────────────────────────────────────────────────────────────

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

test_that("print.taxodist_path runs without error and returns invisibly", {
  mockery::stub(taxo_path, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  })
  result <- taxo_path("Tyrannosaurus", "Triceratops")
  expect_no_error(print(result))
  expect_invisible(print(result))
})

# ── Mock tests ────────────────────────────────────────────────────────────────

test_that("get_taxonomicon_id returns NULL on network failure", {
  clear_cache()
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) NULL)
  expect_warning(result <- get_taxonomicon_id("Tyrannosaurus"), "Cannot reach")
  expect_null(result)
})

test_that("get_taxonomicon_id returns NULL on bad status", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 404L)
  expect_warning(result <- get_taxonomicon_id("Tyrannosaurus"), "Cannot reach")
  expect_null(result)
})

test_that("get_taxonomicon_id returns NULL on an unparseable response", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) "not html")
  mockery::stub(
    get_taxonomicon_id,
    "rvest::read_html",
    function(...) stop("parse failure")
  )
  expect_warning(
    result <- get_taxonomicon_id("Tyrannosaurus"),
    "Could not parse"
  )
  expect_null(result)
})

test_that("get_lineage_by_id returns NULL on network failure", {
  clear_cache()
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) NULL)
  result <- get_lineage_by_id("12345")
  expect_null(result)
})

test_that("get_lineage returns NULL when id not found", {
  clear_cache()
  mockery::stub(get_lineage, "get_taxonomicon_id", function(...) NULL)
  result <- get_lineage("Fakeosaurus")
  expect_null(result)
})

test_that("cache is used on second call to get_taxonomicon_id", {
  clear_cache()
  assign("id_Tyrannosaurus", "50841", envir = taxodist:::.taxodist_cache)
  call_count <- 0L
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) {
    call_count <<- call_count + 1L
    NULL
  })
  result <- get_taxonomicon_id("Tyrannosaurus")
  expect_equal(call_count, 0L)
  expect_equal(result, "50841")
})

test_that("taxo_distance works with mocked lineages", {
  mockery::stub(taxo_distance, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  })
  result <- taxo_distance("Tyrannosaurus", "Triceratops")
  expect_s3_class(result, "taxodist_result")
  expect_equal(result$mrca, "Dinosauria")
})

test_that("closest_relative works with mocked lineages", {
  mockery::stub(closest_relative, "get_lineage", function(taxon, ...) {
    lins <- list(
      Tyrannosaurus = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"),
      Velociraptor  = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Velociraptor"),
      Triceratops   = c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
    )
    lins[[taxon]]
  })
  result <- closest_relative("Tyrannosaurus", c("Velociraptor", "Triceratops"))
  expect_equal(nrow(result), 2)
  expect_equal(result$taxon[1], "Velociraptor")
})

test_that("lineage_depth works with mocked lineage", {
  mockery::stub(lineage_depth, "get_lineage",
                function(...) c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"))
  expect_equal(lineage_depth("Tyrannosaurus"), 5L)
})

test_that("shared_clades works with mocked lineages", {
  mockery::stub(shared_clades, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  })
  result <- shared_clades("Tyrannosaurus", "Triceratops")
  expect_equal(result, c("Biota", "Animalia", "Dinosauria"))
})

test_that("is_member works with mocked lineage", {
  mockery::stub(is_member, "get_lineage",
                function(...) c("Biota", "Animalia", "Dinosauria", "Theropoda"))
  expect_true(is_member("Tyrannosaurus", "Dinosauria"))
  expect_false(is_member("Tyrannosaurus", "Mammalia"))
})

test_that("is_member requires a complete clade name", {
  mockery::stub(is_member, "get_lineage",
                function(...) c("Biota", "Animalia", "Dinosauria", "Theropoda"))
  expect_false(is_member("Tyrannosaurus", "Dino"))
  expect_false(is_member("Tyrannosaurus", "Ther"))
  expect_false(is_member("Tyrannosaurus", "Animal"))
})

test_that("is_member matches case-insensitively and ignores outer whitespace", {
  mockery::stub(is_member, "get_lineage",
                function(...) c("Biota", "Animalia", "Dinosauria", "Theropoda"))
  expect_true(is_member("Tyrannosaurus", "dinosauria"))
  expect_true(is_member("Tyrannosaurus", "  THEROPODA  "))
})

test_that("is_member treats regular expression characters literally", {
  mockery::stub(is_member, "get_lineage",
                function(...) c("Biota", "Clade (example)", "Species+"))
  expect_true(is_member("Example", "Clade (example)"))
  expect_true(is_member("Example", "Species+"))
  expect_false(is_member("Example", "Clade ("))
})

test_that("compare_lineages works with mocked lineages", {
  mockery::stub(compare_lineages, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  })
  result <- compare_lineages("Tyrannosaurus", "Triceratops")
  expect_equal(result$mrca_depth, 3L)
})

test_that("distance_matrix works with mocked lineages", {
  mockery::stub(distance_matrix, "get_lineage", function(taxon, ...) {
    lins <- list(
      Tyrannosaurus = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"),
      Velociraptor  = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Velociraptor"),
      Triceratops   = c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
    )
    lins[[taxon]]
  })
  mat <- distance_matrix(c("Tyrannosaurus", "Velociraptor", "Triceratops"),
                         progress = FALSE)
  expected <- matrix(c(
    0, 1 / 4, 1 / 3,
    1 / 4, 0, 1 / 3,
    1 / 3, 1 / 3, 0
  ), nrow = 3L, byrow = TRUE)
  dimnames(expected) <- rep(list(
    c("Tyrannosaurus", "Velociraptor", "Triceratops")
  ), 2L)

  expect_s3_class(mat, "dist")
  expect_equal(as.matrix(mat), expected)
  expect_equal(attr(mat, "Labels"), colnames(expected))
  expect_equal(attr(mat, "Size"), 3L)
  expect_false(attr(mat, "Diag"))
  expect_false(attr(mat, "Upper"))
})

test_that("distance_matrix returns an empty dist for no taxa", {
  result <- distance_matrix(character(0), progress = FALSE)
  expect_s3_class(result, "dist")
  expect_length(result, 0L)
  expect_equal(attr(result, "Size"), 0L)
})

test_that("distance_matrix handles one taxon without pairwise work", {
  mockery::stub(
    distance_matrix,
    "get_lineage",
    function(...) c("Biota", "Animalia", "Homo")
  )
  result <- distance_matrix("Homo", progress = FALSE)
  expect_s3_class(result, "dist")
  expect_length(result, 0L)
  expect_equal(attr(result, "Size"), 1L)
  expect_equal(attr(result, "Labels"), "Homo")
  expect_equal(as.matrix(result), matrix(
    0,
    nrow = 1L,
    dimnames = list("Homo", "Homo")
  ))
})

test_that("distance_matrix reuses final resolved lineages", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  calls <- 0L
  mockery::stub(distance_matrix, "get_lineage", function(taxon, ...) {
    calls <<- calls + 1L
    c("Biota", "Animalia", taxon)
  })

  taxa <- c("Alpha", "Beta", "Gamma")
  first <- distance_matrix(taxa, progress = FALSE)
  second <- distance_matrix(taxa, progress = FALSE)

  expect_equal(calls, length(taxa))
  expect_equal(second, first)
})

test_that("check_coverage returns named logical vector", {
  mockery::stub(check_coverage, "get_taxonomicon_id",
                function(taxon, ...) if (taxon == "Fakeosaurus") NULL else "12345")
  result <- check_coverage(c("Tyrannosaurus", "Fakeosaurus"))
  expect_type(result, "logical")
  expect_true(result["Tyrannosaurus"])
  expect_false(result["Fakeosaurus"])
})

test_that("taxo_resolve preserves resolved, ambiguous, and unresolved states", {
  search_results <- list(
    Alpha = data.frame(id = "101", name = "Alpha one"),
    Nereis = data.frame(
      id = c("201", "202"),
      name = c("Nereis animal one", "Nereis animal two")
    )
  )
  lineages <- list(
    `101` = c("Biota", "Animalia", "Alpha"),
    `201` = c("Biota", "Animalia", "Nereis"),
    `202` = c("Biota", "Animalia", "Nereis")
  )
  mockery::stub(taxo_resolve, ".taxo_search_details", function(taxon, ...) {
    results <- search_results[[taxon]]
    if (is.null(results)) {
      list(status = "not_found", results = NULL)
    } else {
      list(status = "ok", results = results)
    }
  })
  mockery::stub(taxo_resolve, "get_lineage_by_id", function(taxon_id, ...) {
    lineages[[taxon_id]]
  })

  expect_warning(
    result <- taxo_resolve(
      c("Alpha", "Nereis", "Missing"),
      progress = FALSE
    ),
    "Ambiguous taxon names"
  )

  expect_s3_class(result, "taxodist_resolution")
  expect_equal(result$status, c("resolved", "ambiguous", "unresolved"))
  expect_equal(result$id, c("101", "201", NA_character_))
  expect_equal(result$n_candidates, c(1L, 2L, 0L))
  expect_equal(result$lineage_depth, c(3L, 3L, NA_integer_))
  expect_equal(nrow(result$candidates[[2]]), 2L)
  expect_null(result$lineage[[3]])
})

test_that("taxo_resolve supports explicit ambiguity policies and numeric IDs", {
  candidates <- data.frame(
    id = c("201", "202"),
    name = c("Nereis animal one", "Nereis animal two")
  )
  mockery::stub(taxo_resolve, ".taxo_search_details", function(...) {
    list(status = "ok", results = candidates)
  })
  mockery::stub(taxo_resolve, "get_lineage_by_id", function(taxon_id, ...) {
    c("Biota", "Animalia", if (taxon_id == "999") "Direct" else "Nereis")
  })

  expect_no_warning(
    first <- taxo_resolve("Nereis", ambiguity = "first", progress = FALSE)
  )
  expect_equal(first$status, "ambiguous")
  expect_equal(first$id, "201")

  expect_error(
    taxo_resolve("Nereis", ambiguity = "error", progress = FALSE),
    "Ambiguous taxon names"
  )

  direct <- taxo_resolve("999", progress = FALSE)
  expect_equal(direct$status, "resolved")
  expect_equal(direct$id, "999")
  expect_equal(direct$resolved_name, "Direct")
})

test_that("taxo_resolve validates inputs", {
  expect_error(taxo_resolve(1:3, progress = FALSE), "character vector")
  expect_error(taxo_resolve(c("Alpha", NA), progress = FALSE), "missing or empty")
  expect_error(taxo_resolve(c("Alpha", ""), progress = FALSE), "missing or empty")
  expect_error(
    taxo_resolve("Alpha", ambiguity = "guess", progress = FALSE),
    "arg"
  )
})

test_that("taxo_resolve updates an explicit progress bar", {
  mockery::stub(taxo_resolve, ".taxo_search_details", function(taxon, ...) {
    list(
      status = "ok",
      results = data.frame(id = "101", name = taxon)
    )
  })
  mockery::stub(taxo_resolve, "get_lineage_by_id", function(taxon_id, ...) {
    c("Biota", "Animalia", "Alpha")
  })

  expect_no_error(
    result <- taxo_resolve("Alpha", progress = TRUE)
  )
  expect_equal(result$status, "resolved")
})

test_that("taxo_resolve retrieves duplicate inputs only once", {
  search_calls <- 0L
  lineage_calls <- 0L
  mockery::stub(taxo_resolve, ".taxo_search_details", function(taxon, ...) {
    search_calls <<- search_calls + 1L
    list(
      status = "ok",
      results = data.frame(
        id = if (taxon == "Alpha") "101" else "102",
        name = taxon
      )
    )
  })
  mockery::stub(taxo_resolve, "get_lineage_by_id", function(taxon_id, ...) {
    lineage_calls <<- lineage_calls + 1L
    c("Biota", if (taxon_id == "101") "Alpha" else "Beta")
  })

  result <- taxo_resolve(
    c("Alpha", "Beta", "Alpha", "Beta"),
    progress = FALSE
  )
  expect_equal(search_calls, 2L)
  expect_equal(lineage_calls, 2L)
  expect_equal(result$input, c("Alpha", "Beta", "Alpha", "Beta"))
  expect_equal(result$id, c("101", "102", "101", "102"))
})

test_that("taxo_resolve distinguishes retrieval errors from absent names", {
  mockery::stub(taxo_resolve, ".taxo_search_details", function(taxon, ...) {
    if (taxon == "Offline") {
      list(status = "retrieval_error", results = NULL)
    } else {
      list(status = "not_found", results = NULL)
    }
  })

  result <- taxo_resolve(c("Offline", "Missing"), progress = FALSE)
  expect_equal(result$status, c("retrieval_error", "unresolved"))
  expect_true(all(is.na(result$id)))
  expect_true(all(vapply(result$lineage, is.null, logical(1))))
})

test_that("taxo_resolve preserves candidates when every lineage retrieval fails", {
  candidates <- data.frame(
    id = c("201", "202"),
    name = c("Candidate one", "Candidate two"),
    stringsAsFactors = FALSE
  )
  mockery::stub(taxo_resolve, ".taxo_search_details", function(...) {
    list(status = "ok", results = candidates)
  })
  mockery::stub(taxo_resolve, "get_lineage_by_id", function(...) NULL)

  result <- taxo_resolve("Unavailable", progress = FALSE)
  expect_equal(result$status, "retrieval_error")
  expect_equal(result$n_candidates, 2L)
  expect_equal(result$candidates[[1]], candidates)
  expect_null(result$lineage[[1]])

  bundle <- taxo_bundle(result, progress = FALSE)
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  write_taxodist_bundle(bundle, file)
  restored <- read_taxodist_bundle(file)
  expect_equal(restored$resolution$status, "retrieval_error")
  expect_equal(restored$resolution$n_candidates, 2L)
})

test_that("taxo_resolve marks an unavailable numeric ID as retrieval error", {
  mockery::stub(taxo_resolve, "get_lineage_by_id", function(...) NULL)
  result <- taxo_resolve("999", progress = FALSE)

  expect_equal(result$status, "retrieval_error")
  expect_equal(result$n_candidates, 0L)
  expect_true(is.na(result$id))
  expect_null(result$lineage[[1]])
})

test_that("taxo_from_lineages creates an offline resolution", {
  lineages <- list(
    Alpha = c("Biota", "Animalia", "Alpha"),
    Beta = c("Biota", "Animalia", "Beta")
  )
  resolution <- taxo_from_lineages(lineages, source = "Curated study")

  expect_s3_class(resolution, "taxodist_resolution")
  expect_equal(resolution$input, c("Alpha", "Beta"))
  expect_equal(resolution$resolved_name, c("Alpha", "Beta"))
  expect_equal(resolution$id, c("custom:Alpha", "custom:Beta"))
  expect_equal(resolution$status, c("resolved", "resolved"))
  expect_equal(resolution$lineage_depth, c(3L, 3L))
  expect_equal(attr(resolution, "source"), "Curated study")
  expect_true(is.na(attr(resolution, "source_url")))
  expect_equal(as.matrix(distance_matrix(resolution))[1, 2], 0.5)
})

test_that("taxo_from_lineages accepts named IDs and validates inputs", {
  lineages <- list(
    Alpha = c("Root", "Alpha"),
    Beta = c("Root", "Beta")
  )
  resolution <- taxo_from_lineages(
    lineages,
    ids = c(Beta = "B", Alpha = "A")
  )
  expect_equal(resolution$id, c("A", "B"))

  expect_error(taxo_from_lineages(unname(lineages)), "unique, non-empty names")
  expect_error(
    taxo_from_lineages(list(Alpha = c("Root", NA_character_))),
    "Every lineage"
  )
  expect_error(
    taxo_from_lineages(lineages, ids = c("same", "same")),
    "unique, non-empty"
  )
  expect_error(taxo_from_lineages("not a list"), "named list")
  expect_error(taxo_from_lineages(lineages, source = ""), "source")
  expect_error(taxo_from_lineages(lineages, ids = 1:2), "character vector")
  expect_error(
    taxo_from_lineages(lineages, ids = c(Alpha = "A")),
    "every lineage name"
  )
  expect_error(
    taxo_from_lineages(lineages, ids = "only-one"),
    "one unique"
  )
})

test_that("taxo_bundle preserves custom lineage provenance", {
  resolution <- taxo_from_lineages(list(
    Alpha = c("Biota", "Alpha"),
    Beta = c("Biota", "Beta")
  ), source = "Local taxonomy")
  bundle <- taxo_bundle(resolution, progress = FALSE)

  expect_equal(bundle$source$name, "Local taxonomy")
  expect_true(is.na(bundle$source$url))
  expect_equal(as.matrix(bundle$matrix)[1, 2], 1)

  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  write_taxodist_bundle(bundle, file)
  restored <- read_taxodist_bundle(file)
  expect_equal(restored$source$name, "Local taxonomy")
  expect_true(is.na(restored$source$url))
  expect_equal(attr(restored$resolution, "source"), "Local taxonomy")
  expect_true(is.na(attr(restored$resolution, "source_url")))

  attr(resolution, "source_url") <- NULL
  without_url <- taxo_bundle(resolution, progress = FALSE)
  expect_true(is.na(without_url$source$url))
})

test_that("distance_matrix reuses a taxodist_resolution without retrieval", {
  resolution <- data.frame(
    input = c("Alpha query", "Beta query", "Missing query"),
    resolved_name = c("Alpha", "Beta", NA_character_),
    id = c("101", "102", NA_character_),
    status = c("resolved", "resolved", "unresolved"),
    n_candidates = c(1L, 1L, 0L),
    lineage_depth = c(3L, 3L, NA_integer_),
    stringsAsFactors = FALSE
  )
  resolution$lineage <- I(list(
    c("Biota", "Animalia", "Alpha"),
    c("Biota", "Animalia", "Beta"),
    NULL
  ))
  resolution$candidates <- I(list(
    data.frame(id = "101", name = "Alpha"),
    data.frame(id = "102", name = "Beta"),
    data.frame(id = character(), name = character())
  ))
  class(resolution) <- c("taxodist_resolution", "data.frame")

  matrix <- as.matrix(distance_matrix(resolution, progress = FALSE))
  expect_equal(rownames(matrix), resolution$input)
  expect_equal(colnames(matrix), resolution$input)
  expect_equal(matrix[1, 2], 0.5)
  expect_true(is.na(matrix[1, 3]))
  expect_true(is.na(matrix[2, 3]))
})

test_that("taxo_bundle combines resolution, matrix, metric, and provenance", {
  resolution <- data.frame(
    input = c("Alpha query", "Beta query", "Missing query"),
    resolved_name = c("Alpha", "Beta", NA_character_),
    id = c("101", "102", NA_character_),
    status = c("resolved", "ambiguous", "unresolved"),
    n_candidates = c(1L, 2L, 0L),
    lineage_depth = c(3L, 3L, NA_integer_),
    stringsAsFactors = FALSE
  )
  resolution$lineage <- I(list(
    c("Biota", "Animalia", "Alpha"),
    c("Biota", "Animalia", "Beta"),
    NULL
  ))
  resolution$candidates <- I(list(
    data.frame(id = "101", name = "Alpha"),
    data.frame(id = c("102", "103"), name = c("Beta one", "Beta two")),
    data.frame(id = character(), name = character())
  ))
  class(resolution) <- c("taxodist_resolution", "data.frame")
  attr(resolution, "source") <- "The Taxonomicon"
  attr(resolution, "retrieved_at") <- "2026-09-23 12:00:00 UTC"

  bundle <- taxo_bundle(resolution, progress = FALSE)
  expect_s3_class(bundle, "taxodist_bundle")
  expect_equal(bundle$schema_version, "1.0")
  expect_equal(bundle$source$name, "The Taxonomicon")
  expect_equal(bundle$metric$name, "inverse_mrca_depth")
  expect_identical(bundle$resolution, resolution)
  expect_s3_class(bundle$matrix, "dist")
  expect_equal(attr(bundle$matrix, "Labels"), resolution$input)
  expect_identical(distance_matrix(bundle), bundle$matrix)

  invalid <- bundle
  invalid$resolution$status[[1]] <- "unresolved"
  expect_error(
    write_taxodist_bundle(invalid, tempfile(fileext = ".json")),
    "incomplete unresolved record"
  )
})

test_that("taxo_bundle supports an empty character input and JSON round-trip", {
  bundle <- taxo_bundle(character(0), progress = FALSE)
  expect_s3_class(bundle, "taxodist_bundle")
  expect_equal(nrow(bundle$resolution), 0L)
  expect_equal(attr(bundle$matrix, "Size"), 0L)

  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  write_taxodist_bundle(bundle, file)
  restored <- read_taxodist_bundle(file)
  expect_equal(nrow(restored$resolution), 0L)
  expect_equal(attr(restored$matrix, "Size"), 0L)
})

test_that("taxo_bundle supplies fallback provenance and package metadata", {
  resolution <- taxo_from_lineages(list(
    Alpha = c("Biota", "Alpha")
  ))
  attr(resolution, "source") <- ""
  attr(resolution, "source_url") <- NULL
  mockery::stub(taxo_bundle, "utils::packageVersion", function(...) {
    stop("version unavailable")
  })

  bundle <- taxo_bundle(resolution, progress = FALSE)
  expect_true(is.na(bundle$software$version))
  expect_equal(bundle$source$name, "The Taxonomicon")
  expect_equal(bundle$source$url, "http://taxonomicon.taxonomy.nl")
})

test_that("bundle validation rejects each malformed component", {
  base <- taxo_bundle(taxo_from_lineages(list(
    Alpha = c("Biota", "Animalia", "Alpha"),
    Beta = c("Biota", "Animalia", "Beta")
  )), progress = FALSE)

  expect_bad <- function(change, pattern) {
    invalid <- change(base)
    expect_error(
      taxodist:::.validate_taxodist_bundle(invalid),
      pattern
    )
  }

  expect_bad(function(x) {
    x$metric <- NULL
    x
  }, "required fields")
  expect_bad(function(x) {
    x$schema_version <- "2.0"
    x
  }, "Unsupported")
  expect_bad(function(x) {
    class(x$resolution) <- "data.frame"
    x
  }, "wrong class")
  expect_bad(function(x) {
    x$resolution$candidates <- NULL
    x
  }, "resolution fields")
  expect_bad(function(x) {
    x$resolution$status[[1]] <- "mystery"
    x
  }, "unknown resolution status")
  expect_bad(function(x) {
    x$resolution$n_candidates[[1]] <- NA_integer_
    x
  }, "candidate counts")
  expect_bad(function(x) {
    x$resolution$candidates[1] <- list(list(id = "x", name = "x"))
    x
  }, "malformed candidate table")
  expect_bad(function(x) {
    x$resolution$n_candidates[[1]] <- 2L
    x
  }, "candidate count mismatch")
  expect_bad(function(x) {
    x$resolution$lineage[1] <- list(NULL)
    x
  }, "lineage depth mismatch")
  expect_bad(function(x) {
    x$resolution$lineage[1] <- list(1:3)
    x
  }, "malformed lineage")
  expect_bad(function(x) {
    x$resolution$lineage_depth[[1]] <- 99L
    x
  }, "lineage depth mismatch")
  expect_bad(function(x) {
    x$resolution$id[[1]] <- NA_character_
    x
  }, "incomplete resolved record")
  expect_bad(function(x) {
    x$resolution$candidates[[1]]$id[[1]] <- "different"
    x
  }, "selected candidate mismatch")
  expect_bad(function(x) {
    x$resolution$status[[1]] <- "unresolved"
    x
  }, "incomplete unresolved record")
  expect_bad(function(x) {
    x$resolution$status[[1]] <- "unresolved"
    x$resolution$id[[1]] <- NA_character_
    x$resolution$resolved_name[[1]] <- NA_character_
    x$resolution$lineage[1] <- list(NULL)
    x$resolution$lineage_depth[[1]] <- NA_integer_
    x
  }, "unresolved record has candidates")
  expect_bad(function(x) {
    class(x$matrix) <- NULL
    x
  }, "must be a.*dist")
  expect_bad(function(x) {
    attr(x$matrix, "Labels") <- rev(attr(x$matrix, "Labels"))
    x
  }, "matrix labels")
  expect_bad(function(x) {
    x$matrix[] <- 0.75
    x
  }, "stored distances")
})

test_that("taxodist bundle JSON round-trip preserves scientific contents", {
  resolution <- data.frame(
    input = c("Alpha query", "Missing query"),
    resolved_name = c("Alpha", NA_character_),
    id = c("101", NA_character_),
    status = c("resolved", "unresolved"),
    n_candidates = c(1L, 0L),
    lineage_depth = c(3L, NA_integer_),
    stringsAsFactors = FALSE
  )
  resolution$lineage <- I(list(
    c("Biota", "Animalia", "Alpha"),
    NULL
  ))
  resolution$candidates <- I(list(
    data.frame(id = "101", name = "Alpha accepted"),
    data.frame(id = character(), name = character())
  ))
  class(resolution) <- c("taxodist_resolution", "data.frame")
  attr(resolution, "source") <- "The Taxonomicon"
  attr(resolution, "retrieved_at") <- "2026-09-23 12:00:00 UTC"
  bundle <- taxo_bundle(resolution, progress = FALSE)
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)

  expect_invisible(path <- write_taxodist_bundle(bundle, file))
  expect_true(file.exists(file))
  expect_true(jsonlite::validate(paste(readLines(file, warn = FALSE), collapse = "\n")))

  restored <- read_taxodist_bundle(file)
  expect_s3_class(restored, "taxodist_bundle")
  expect_equal(restored$schema_version, bundle$schema_version)
  expect_equal(restored$resolution$input, resolution$input)
  expect_equal(restored$resolution$status, resolution$status)
  expect_equal(restored$resolution$id, resolution$id)
  expect_equal(restored$resolution$lineage, resolution$lineage)
  expect_equal(restored$resolution$candidates, resolution$candidates)
  expect_equal(as.matrix(restored$matrix), as.matrix(bundle$matrix))
  expect_equal(path, normalizePath(file, winslash = "/"))
})

test_that("bundle JSON preserves one-element arrays and infinite distances", {
  resolution <- data.frame(
    input = "Disconnected",
    resolved_name = "Disconnected",
    id = "500",
    status = "resolved",
    n_candidates = 1L,
    lineage_depth = 1L,
    stringsAsFactors = FALSE
  )
  resolution$lineage <- I(list("Disconnected"))
  resolution$candidates <- I(list(
    data.frame(id = "500", name = "Disconnected")
  ))
  class(resolution) <- c("taxodist_resolution", "data.frame")
  bundle <- taxo_bundle(resolution, progress = FALSE)
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  write_taxodist_bundle(bundle, file)

  raw <- jsonlite::read_json(file, simplifyVector = FALSE)
  expect_type(raw$matrix$labels, "list")
  expect_length(raw$matrix$labels, 1L)
  expect_type(raw$taxa[[1]]$lineage, "list")
  expect_length(raw$taxa[[1]]$lineage, 1L)

  two <- resolution[c(1, 1), , drop = FALSE]
  two$input <- c("Disconnected A", "Disconnected B")
  two$id <- c("500", "501")
  two$lineage <- I(list("Root A", "Root B"))
  two$candidates <- I(list(
    data.frame(id = "500", name = "Disconnected A"),
    data.frame(id = "501", name = "Disconnected B")
  ))
  class(two) <- c("taxodist_resolution", "data.frame")
  infinite_bundle <- taxo_bundle(two, progress = FALSE)
  expect_true(is.infinite(as.matrix(infinite_bundle$matrix)[1, 2]))
  write_taxodist_bundle(infinite_bundle, file)
  restored <- read_taxodist_bundle(file)
  expect_equal(as.matrix(restored$matrix)[1, 2], Inf)
})

test_that("bundle readers and writers reject invalid inputs", {
  expect_error(
    write_taxodist_bundle(list(), tempfile(fileext = ".json")),
    "taxodist_bundle"
  )
  expect_error(
    read_taxodist_bundle("missing-taxodist-bundle.json"),
    "not found"
  )

  bad <- tempfile(fileext = ".json")
  on.exit(unlink(bad), add = TRUE)
  writeLines('{"format":"something_else","schema_version":"1.0"}', bad)
  expect_error(read_taxodist_bundle(bad), "unrecognized format")

  writeLines("{", bad)
  expect_error(read_taxodist_bundle(bad), "Could not parse")

  writeLines(
    '{"format":"taxodist_bundle","schema_version":"2.0"}',
    bad
  )
  expect_error(read_taxodist_bundle(bad), "Unsupported")

  writeLines(
    '{"format":"taxodist_bundle","schema_version":"1.0"}',
    bad
  )
  expect_error(read_taxodist_bundle(bad), "required fields")
})

test_that("bundle matrix JSON validation covers malformed matrices", {
  matrix_from_json <- taxodist:::.matrix_from_json

  expect_error(
    matrix_from_json(list(
      labels = list("A", "B"),
      values = list(list(0, 1))
    )),
    "row count"
  )
  expect_error(
    matrix_from_json(list(
      labels = list("A", "B"),
      values = list(list(0), list(1, 0))
    )),
    "square"
  )
  expect_error(
    matrix_from_json(list(
      labels = list("A", "B"),
      values = list(list(1, 0.5), list(0.5, 0))
    )),
    "diagonal"
  )
  expect_error(
    matrix_from_json(list(
      labels = list("A", "B"),
      values = list(list(0, 0.25), list(0.5, 0))
    )),
    "symmetric"
  )

  negative_infinity <- matrix_from_json(list(
    labels = list("A", "B"),
    values = list(
      list(0, "-Infinity"),
      list("-Infinity", 0)
    )
  ))
  expect_equal(as.matrix(negative_infinity)[1, 2], -Inf)
})

test_that("print.taxodist_bundle returns invisibly", {
  resolution <- data.frame(
    input = "Alpha",
    resolved_name = "Alpha",
    id = "101",
    status = "resolved",
    n_candidates = 1L,
    lineage_depth = 2L,
    stringsAsFactors = FALSE
  )
  resolution$lineage <- I(list(c("Biota", "Alpha")))
  resolution$candidates <- I(list(data.frame(id = "101", name = "Alpha")))
  class(resolution) <- c("taxodist_resolution", "data.frame")
  bundle <- taxo_bundle(resolution, progress = FALSE)

  expect_no_error(print(bundle))
  expect_invisible(print(bundle))
})

test_that("taxo_distance returns NULL when taxon_a not found", {
  mockery::stub(taxo_distance, "get_lineage", function(taxon, ...) NULL)
  result <- taxo_distance("Fakeosaurus", "Carnotaurus")
  expect_null(result)
})

test_that("taxo_distance returns NULL when taxon_b not found", {
  mockery::stub(taxo_distance, "get_lineage", function(taxon, ...) {
    if (taxon == "Carnotaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Abelisauridae", "Carnotaurus")
    else
      NULL
  })
  result <- taxo_distance("Carnotaurus", "Fakeosaurus")
  expect_null(result)
})

test_that("mrca returns NULL when taxo_distance fails", {
  mockery::stub(mrca, "taxo_distance", function(...) NULL)
  result <- mrca("Fakeosaurus", "Carnotaurus")
  expect_null(result)
})

test_that("mrca returns correct value when taxo_distance succeeds", {
  mockery::stub(mrca, "taxo_distance", function(...) list(mrca = "Dinosauria"))
  result <- mrca("Carnotaurus", "Triceratops")
  expect_equal(result, "Dinosauria")
})

test_that("closest_relative returns NULL when query lineage not found", {
  mockery::stub(closest_relative, "get_lineage", function(...) NULL)
  result <- closest_relative("Fakeosaurus", c("Carnotaurus", "Velociraptor"))
  expect_null(result)
})

test_that("closest_relative returns an empty data frame for no candidates", {
  mockery::stub(closest_relative, "get_lineage", function(...) {
    c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus")
  })
  result <- closest_relative("Tyrannosaurus", character(0))
  expect_s3_class(result, "data.frame")
  expect_named(result, c("taxon", "distance"))
  expect_equal(nrow(result), 0L)
})

test_that("closest_relative handles NULL candidate lineage", {
  mockery::stub(closest_relative, "get_lineage", function(taxon, ...) {
    if (taxon == "Carnotaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Abelisauridae", "Carnotaurus")
    else if (taxon == "Velociraptor")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Dromaeosauridae", "Velociraptor")
    else
      NULL
  })
  result <- closest_relative("Carnotaurus", c("Velociraptor", "Fakeosaurus"))
  expect_equal(nrow(result), 2)
  expect_true(is.na(result$distance[result$taxon == "Fakeosaurus"]))
})

test_that("distance_matrix handles NULL lineage for one taxon", {
  mockery::stub(distance_matrix, "get_lineage", function(taxon, ...) {
    if (taxon == "Fakeosaurus") NULL
    else if (taxon == "Carnotaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Abelisauridae", "Carnotaurus")
    else
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Dromaeosauridae", "Velociraptor")
  })
  mat <- distance_matrix(c("Carnotaurus", "Velociraptor", "Fakeosaurus"),
                         progress = FALSE)
  m <- as.matrix(mat)
  expect_true(is.na(m["Carnotaurus", "Fakeosaurus"]))
  expect_false(is.na(m["Carnotaurus", "Velociraptor"]))
})

test_that("distance_matrix with progress = TRUE runs without error", {
  mockery::stub(distance_matrix, "get_lineage", function(taxon, ...) {
    lins <- list(
      Carnotaurus  = c("Biota", "Animalia", "Dinosauria", "Theropoda",
                       "Abelisauridae", "Carnotaurus"),
      Velociraptor = c("Biota", "Animalia", "Dinosauria", "Theropoda",
                       "Dromaeosauridae", "Velociraptor"),
      Triceratops  = c("Biota", "Animalia", "Dinosauria",
                       "Ornithischia", "Triceratops")
    )
    lins[[taxon]]
  })
  expect_no_error(
    distance_matrix(c("Carnotaurus", "Velociraptor", "Triceratops"),
                    progress = TRUE)
  )
})

test_that("get_taxonomicon_id verbose prints messages on cache hit", {
  clear_cache()
  assign("id_Carnotaurus", "99999", envir = taxodist:::.taxodist_cache)
  expect_no_error(get_taxonomicon_id("Carnotaurus", verbose = TRUE))
})

test_that("get_taxonomicon_id verbose prints warning on network failure", {
  clear_cache()
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) NULL)
  expect_warning(get_taxonomicon_id("Drosophila", verbose = TRUE), "Cannot reach")
})

test_that("get_taxonomicon_id verbose prints warning on bad status", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 503L)
  expect_warning(get_taxonomicon_id("Drosophila", verbose = TRUE), "Cannot reach")
})

test_that("get_lineage_by_id verbose prints messages on cache hit", {
  clear_cache()
  assign("lin_99999",
         c("Biota", "Animalia", "Dinosauria", "Abelisauridae", "Carnotaurus"),
         envir = taxodist:::.taxodist_cache)
  expect_no_error(get_lineage_by_id("99999", verbose = TRUE))
})

test_that("get_lineage_by_id verbose prints warning on network failure", {
  clear_cache()
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) NULL)
  expect_no_error(get_lineage_by_id("00000", verbose = TRUE))
})

test_that("get_lineage_by_id returns NULL when lineage is empty after cleaning", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content",
                function(...) "<html><body></body></html>")
  mockery::stub(get_lineage_by_id, "rvest::read_html",
                function(...) xml2::read_html("<html><body></body></html>"))
  result <- get_lineage_by_id("empty_page")
  expect_null(result)
})

test_that("get_lineage_by_id returns NULL when the response cannot be parsed", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) "not html")
  mockery::stub(
    get_lineage_by_id,
    "rvest::read_html",
    function(...) stop("parse failure")
  )
  expect_null(get_lineage_by_id("12345", verbose = TRUE))
})

test_that("get_lineage_by_id returns NULL on bad HTTP status", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 404L)
  result <- get_lineage_by_id("99999")
  expect_null(result)
})

test_that("get_lineage_by_id verbose prints warning on bad status", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 503L)
  expect_no_error(get_lineage_by_id("99999", verbose = TRUE))
})

test_that("get_lineage_by_id cache hit with verbose prints message", {
  clear_cache()
  assign("lin_50841",
         c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  result <- get_lineage_by_id("50841", verbose = TRUE)
  expect_type(result, "character")
})

test_that("get_lineage_by_id with clean = FALSE returns NULL on network failure", {
  clear_cache()
  assign("lin_clean_test", NULL, envir = taxodist:::.taxodist_cache)
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) NULL)
  result <- get_lineage_by_id("clean_test", clean = FALSE)
  expect_null(result)
})

test_that("print.taxodist_result displays output correctly", {
  lin_a <- c("Biota", "Animalia", "Dinosauria", "Theropoda",
             "Abelisauridae", "Carnotaurus")
  lin_b <- c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
  result <- taxodist:::.compute_distance(lin_a, lin_b, "Carnotaurus", "Triceratops")
  expect_no_error(print(result))
  expect_invisible(print(result))
})

test_that("compare_lineages returns invisible NULL when lineage missing", {
  mockery::stub(compare_lineages, "get_lineage", function(...) NULL)
  result <- compare_lineages("Fakeosaurus", "Carnotaurus")
  expect_null(result)
})

test_that("compare_lineages handles mrca_depth == 0", {
  mockery::stub(compare_lineages, "get_lineage", function(taxon, ...) {
    if (taxon == "Drosophila")
      c("Biota", "Animalia", "Arthropoda", "Insecta", "Drosophila")
    else
      c("Fungi", "Ascomycota", "Saccharomyces")
  })
  result <- compare_lineages("Drosophila", "Saccharomyces")
  expect_equal(result$mrca_depth, 0L)
})

test_that("compare_lineages handles case where one lineage is subset of other", {
  mockery::stub(compare_lineages, "get_lineage", function(taxon, ...) {
    if (taxon == "Dinosauria")
      c("Biota", "Animalia", "Dinosauria")
    else
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Carnotaurus")
  })
  result <- compare_lineages("Dinosauria", "Carnotaurus")
  expect_equal(result$mrca_depth, 3L)
})

test_that("shared_clades returns character(0) when no common ancestor", {
  mockery::stub(shared_clades, "get_lineage", function(taxon, ...) {
    if (taxon == "Drosophila")
      c("Biota", "Animalia", "Arthropoda", "Insecta", "Drosophila")
    else
      c("Fungi", "Ascomycota", "Saccharomyces")
  })
  result <- shared_clades("Drosophila", "Saccharomyces")
  expect_equal(result, character(0))
})

test_that("shared_clades returns NULL when one lineage missing", {
  mockery::stub(shared_clades, "get_lineage", function(...) NULL)
  result <- shared_clades("Fakeosaurus", "Carnotaurus")
  expect_null(result)
})

test_that("is_member returns NULL when lineage not found", {
  mockery::stub(is_member, "get_lineage", function(...) NULL)
  result <- is_member("Fakeosaurus", "Dinosauria")
  expect_null(result)
})

test_that("filter_clade handles NULL result from is_member", {
  mockery::stub(filter_clade, "is_member", function(taxon, clade, ...) {
    if (taxon == "Fakeosaurus") NULL
    else clade %in% list(
      Carnotaurus = c("Dinosauria", "Theropoda"),
      Drosophila  = c("Animalia", "Insecta")
    )[[taxon]]
  })
  result <- filter_clade(
    c("Carnotaurus", "Fakeosaurus", "Drosophila"), "Dinosauria"
  )
  expect_equal(result, "Carnotaurus")
})

test_that("lineage_depth returns NULL when lineage not found", {
  mockery::stub(lineage_depth, "get_lineage", function(...) NULL)
  result <- lineage_depth("Fakeosaurus")
  expect_null(result)
})

test_that("get_taxonomicon_id parses HTML and returns id", {
  clear_cache()
  fake_html <- '
    <html><body><table>
      <tr>
        <td>Carnotaurus - animal - dinosaur</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=12345&src=0">tree</a></td>
      </tr>
    </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)
  mockery::stub(get_taxonomicon_id, "get_lineage_by_id", function(...) c("Biota", "Animalia"))
  result <- get_taxonomicon_id("Carnotaurus")
  expect_equal(result, "12345")
})

test_that("get_taxonomicon_id skips astronomical entries", {
  clear_cache()
  fake_html <- '
    <html><body><table>
      <tr>
        <td>Carnotaurus - asteroid - Minor planet</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=99999&src=0">tree</a></td>
      </tr>
      <tr>
        <td>Carnotaurus - animal - dinosaur</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=12345&src=0">tree</a></td>
      </tr>
    </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)
  mockery::stub(get_taxonomicon_id, "get_lineage_by_id", function(...) c("Biota", "Animalia"))
  result <- get_taxonomicon_id("Carnotaurus", verbose = TRUE)
  expect_equal(result, "12345")
})

test_that("get_taxonomicon_id skips row matching both bio and astronomical keywords", {
  clear_cache()
  fake_html <- '
    <html><body><table>
      <tr>
        <td>Pterodactylus - animal - Minor planet asteroid</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=99999&src=0">wrong</a></td>
      </tr>
      <tr>
        <td>Pterodactylus - animal - reptile</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=42042&src=0">tree</a></td>
      </tr>
    </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)
  mockery::stub(get_taxonomicon_id, "get_lineage_by_id", function(...) c("Biota", "Animalia"))
  result <- get_taxonomicon_id("Pterodactylus", verbose = TRUE)
  expect_equal(result, "42042")
})

test_that("get_taxonomicon_id returns NULL when bio row has no tree link", {
  clear_cache()
  fake_html <- '
    <html><body><table>
      <tr>
        <td>Quercus - plant</td>
        <td><a href="SomeOtherPage.aspx?id=999">no tree link</a></td>
      </tr>
    </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)
  result <- get_taxonomicon_id("Quercus", verbose = TRUE)
  expect_null(result)
})

test_that("get_lineage_by_id parses HTML and returns lineage", {
  clear_cache()
  fake_html <- '
    <html><body>
      <a href="TaxonTree.aspx?id=1&src=0">Biota</a>
      <a href="TaxonTree.aspx?id=2&src=0">Animalia</a>
      <a href="TaxonTree.aspx?id=3&src=0">Dinosauria</a>
      <a href="TaxonTree.aspx?id=4&src=0">Theropoda</a>
      <a href="TaxonTree.aspx?id=5&src=0">Carnotaurus</a>
    </body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) fake_html)
  result <- get_lineage_by_id("12345", verbose = TRUE)
  expect_type(result, "character")
  expect_true("Dinosauria" %in% result)
  expect_true("Carnotaurus" %in% result)
})

test_that("get_lineage passes clean and verbose through to get_lineage_by_id", {
  clear_cache()
  mockery::stub(get_lineage, "get_taxonomicon_id", function(...) "12345")
  mockery::stub(get_lineage, "get_lineage_by_id",
                function(id, clean, verbose) {
                  expect_equal(id, "12345")
                  expect_false(clean)
                  expect_true(verbose)
                  c("Biota", "Animalia", "Plantae", "Quercus")
                })
  result <- get_lineage("Quercus", clean = FALSE, verbose = TRUE)
  expect_equal(result, c("Biota", "Animalia", "Plantae", "Quercus"))
})

test_that("plot.taxodist_cluster runs without error", {
  m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
              nrow = 3,
              dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  cl <- taxo_cluster(stats::as.dist(m))
  expect_no_error(plot(cl))
})

test_that("plot.taxodist_ord runs without error", {
  m <- matrix(c(0, 0.2, 0.5, 0.2, 0, 0.3, 0.5, 0.3, 0),
              nrow = 3,
              dimnames = list(c("A","B","C"), c("A","B","C")))
  ord <- taxo_ordinate(stats::as.dist(m))
  expect_no_error(plot(ord))
  expect_invisible(plot(ord))
})

# ── Mocks for taxo_search ─────────────────────────────────────────────────────

test_that("user agent has a development fallback", {
  user_agent <- taxodist:::.taxodist_user_agent
  mockery::stub(user_agent, "utils::packageVersion", function(...) {
    stop("version unavailable")
  })
  expect_equal(user_agent(), "taxodist R package development")
})

test_that("taxo_search keeps diagnostics out of the public API", {
  expect_named(formals(taxo_search), c("taxon", "verbose"))
  expected <- data.frame(id = "101", name = "Alpha")
  mockery::stub(taxo_search, ".taxo_search_details", function(...) {
    list(status = "ok", results = expected)
  })
  expect_equal(taxo_search("Alpha"), expected)
})

test_that("taxo_search details report network, status, and parsing failures", {
  clear_cache()
  search_details <- taxodist:::.taxo_search_details
  mockery::stub(search_details, "httr::GET", function(...) stop("Network error"))
  expect_equal(
    search_details("Bacteria", verbose = TRUE)$status,
    "retrieval_error"
  )

  fake_response <- structure(list(), class = "response")
  mockery::stub(search_details, "httr::GET", function(...) fake_response)
  mockery::stub(search_details, "httr::status_code", function(...) 503L)
  expect_equal(search_details("Bacteria")$status, "retrieval_error")

  mockery::stub(search_details, "httr::status_code", function(...) 200L)
  mockery::stub(search_details, "httr::content", function(...) "not html")
  mockery::stub(search_details, "rvest::read_html", function(...) {
    stop("parse failure")
  })
  expect_equal(search_details("Bacteria")$status, "retrieval_error")
})

test_that("taxo_search returns NULL when no matches are found", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  search_details <- taxodist:::.taxo_search_details
  mockery::stub(search_details, "httr::GET", function(...) fake_response)
  mockery::stub(search_details, "httr::status_code", function(...) 200L)
  mockery::stub(
    search_details,
    "httr::content",
    function(...) "<html><body><table></table></body></html>"
  )

  details <- search_details("EmptyTaxon", verbose = TRUE)
  expect_equal(details$status, "not_found")
  expect_null(details$results)
})

test_that("taxo_search parses HTML, applies skips, dedups, and returns data.frame", {
  clear_cache()
  fake_html <- '
    <html><body><table>
      <tr>
        <td>Astronomical planet asteroid</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=111">ignore</a></td>
      </tr>
      <tr>
        <td>No links here</td>
        <td>Just text</td>
      </tr>
      <tr>
        <td>Invalid class</td>
        <td><a class="Invalid" href="TaxonTree.aspx?id=222">ignore</a></td>
      </tr>
      <tr>
        <td>Missing ID</td>
        <td><a class="Valid" href="TaxonTree.aspx?wrong=333">ignore</a></td>
      </tr>
      <tr>
        <td><a class="Valid" href="TaxonTree.aspx?id=444">N|T|P|R|B|L Bacteria (Kingdom)</a></td>
      </tr>
      <tr>
        <td><a class="Valid" href="TaxonTree.aspx?id=444">N|T|P|R|B|L Bacteria (Kingdom) Duplicated</a></td>
      </tr>
      <tr>
        <td><a class="Valid" href="TaxonTree.aspx?id=555">N|T|P|R|B|L Bacteria (Domain)</a></td>
      </tr>
    </table></body></html>'

  fake_response <- structure(list(), class = "response")
  search_details <- taxodist:::.taxo_search_details
  mockery::stub(search_details, "httr::GET", function(...) fake_response)
  mockery::stub(search_details, "httr::status_code", function(...) 200L)
  mockery::stub(search_details, "httr::content", function(...) fake_html)

  details <- search_details("Bacteria", verbose = TRUE)
  df <- details$results

  expect_equal(details$status, "ok")
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 2)
  expect_equal(df$id, c("444", "555"))
  expect_equal(df$name[1], "Bacteria (Kingdom)")
  expect_equal(df$name[2], "Bacteria (Domain)")
})

# ── cache_info ────────────────────────────────────────────────────────────────

test_that("cache_info returns invisible list with correct structure", {
  clear_cache()
  result <- cache_info()
  expect_invisible(cache_info())
  expect_type(result, "list")
  expect_named(result, c("n_lineages", "n_ids", "taxa", "size_bytes"))
})

test_that("cache_info reports zero counts on empty cache", {
  clear_cache()
  result <- cache_info()
  expect_equal(result$n_lineages, 0L)
  expect_equal(result$n_ids, 0L)
  expect_equal(result$taxa, character(0))
})

test_that("cache_info counts lineages and IDs correctly", {
  clear_cache()
  assign("id_Tyrannosaurus", "50841", envir = taxodist:::.taxodist_cache)
  assign("id_Velociraptor",  "12345", envir = taxodist:::.taxodist_cache)
  assign("lin_50841", c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  result <- cache_info()
  expect_equal(result$n_lineages, 1L)
  expect_equal(result$n_ids, 2L)
})

test_that("cache_info taxa names strip the lin_ prefix", {
  clear_cache()
  assign("lin_50841", c("Biota", "Animalia", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  assign("lin_12345", c("Biota", "Animalia", "Velociraptor"),
         envir = taxodist:::.taxodist_cache)
  result <- cache_info()
  expect_setequal(result$taxa, c("50841", "12345"))
})

test_that("cache_info size_bytes is numeric and positive when cache is populated", {
  clear_cache()
  assign("lin_50841", c("Biota", "Animalia", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  result <- cache_info()
  expect_type(result$size_bytes, "double")
  expect_gt(result$size_bytes, 0)
})

test_that("cache_info prints without error", {
  clear_cache()
  assign("lin_50841", c("Biota", "Animalia", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  expect_no_error(cache_info())
})

# ── focal_distances ───────────────────────────────────────────────────────────

test_that("focal_distances returns NULL when focal taxon not found", {
  mockery::stub(focal_distances, "get_lineage", function(...) NULL)
  result <- focal_distances("Fakeosaurus", c("Velociraptor", "Triceratops"))
  expect_null(result)
})

test_that("focal_distances returns a typed empty result for no community", {
  mockery::stub(focal_distances, "get_lineage", function(...) {
    c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus")
  })
  result <- focal_distances(
    "Tyrannosaurus", character(0), progress = FALSE
  )
  expect_s3_class(result, "taxodist_focal")
  expect_s3_class(result, "data.frame")
  expect_named(result, c("taxon", "distance", "mrca", "mrca_depth"))
  expect_equal(nrow(result), 0L)
  expect_equal(attr(result, "focal"), "Tyrannosaurus")
})

test_that("focal_distances returns correct S3 class", {
  mockery::stub(focal_distances, "get_lineage", function(taxon, ...) {
    lins <- list(
      Tyrannosaurus = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"),
      Velociraptor  = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Velociraptor"),
      Triceratops   = c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
    )
    lins[[taxon]]
  })
  result <- focal_distances("Tyrannosaurus",
                            c("Velociraptor", "Triceratops"),
                            progress = FALSE)
  expect_s3_class(result, "taxodist_focal")
  expect_s3_class(result, "data.frame")
})

test_that("focal_distances result has correct columns", {
  mockery::stub(focal_distances, "get_lineage", function(taxon, ...) {
    lins <- list(
      Tyrannosaurus = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"),
      Velociraptor  = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Velociraptor")
    )
    lins[[taxon]]
  })
  result <- focal_distances("Tyrannosaurus", "Velociraptor", progress = FALSE)
  expect_named(result, c("taxon", "distance", "mrca", "mrca_depth"))
})

test_that("focal_distances is sorted by distance ascending", {
  mockery::stub(focal_distances, "get_lineage", function(taxon, ...) {
    lins <- list(
      Tyrannosaurus = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"),
      Velociraptor  = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Velociraptor"),
      Triceratops   = c("Biota", "Animalia", "Dinosauria", "Ornithischia", "Triceratops")
    )
    lins[[taxon]]
  })
  result <- focal_distances("Tyrannosaurus",
                            c("Triceratops", "Velociraptor"),
                            progress = FALSE)
  expect_equal(result$taxon[1], "Velociraptor")
})

test_that("focal_distances handles focal in community with distance 0", {
  mockery::stub(focal_distances, "get_lineage", function(taxon, ...) {
    c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
  })
  result <- focal_distances("Tyrannosaurus",
                            c("Tyrannosaurus", "Tyrannosaurus"),
                            progress = FALSE)
  expect_true(all(result$distance == 0))
  expect_true(all(result$mrca == "Tyrannosaurus"))
})

test_that("focal_distances handles NULL candidate lineage with NA row", {
  mockery::stub(focal_distances, "get_lineage", function(taxon, ...) {
    if (taxon == "Tyrannosaurus")
      c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
    else
      NULL
  })
  result <- focal_distances("Tyrannosaurus", "Fakeosaurus", progress = FALSE)
  expect_true(is.na(result$distance[result$taxon == "Fakeosaurus"]))
  expect_true(is.na(result$mrca[result$taxon == "Fakeosaurus"]))
})

test_that("focal_distances preserves focal attribute", {
  mockery::stub(focal_distances, "get_lineage", function(taxon, ...) {
    c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus")
  })
  result <- focal_distances("Tyrannosaurus", "Tyrannosaurus", progress = FALSE)
  expect_equal(attr(result, "focal"), "Tyrannosaurus")
})

test_that("focal_distances with progress = TRUE runs without error", {
  mockery::stub(focal_distances, "get_lineage", function(taxon, ...) {
    lins <- list(
      Tyrannosaurus = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"),
      Velociraptor  = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Velociraptor")
    )
    lins[[taxon]]
  })
  expect_no_error(
    focal_distances("Tyrannosaurus", "Velociraptor", progress = FALSE)
  )
})

test_that("focal_distances progress bar lines are executed", {
  clear_cache()
  assign("id_Tyrannosaurus", "50841", envir = taxodist:::.taxodist_cache)
  assign("id_Velociraptor",  "12345", envir = taxodist:::.taxodist_cache)
  assign("lin_50841", c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  assign("lin_12345", c("Biota", "Animalia", "Dinosauria", "Theropoda", "Velociraptor"),
         envir = taxodist:::.taxodist_cache)

  expect_no_error(
    focal_distances("Tyrannosaurus", "Velociraptor", progress = TRUE)
  )
})

test_that("print.taxodist_focal runs without error and returns invisibly", {
  mockery::stub(focal_distances, "get_lineage", function(taxon, ...) {
    lins <- list(
      Tyrannosaurus = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"),
      Velociraptor  = c("Biota", "Animalia", "Dinosauria", "Theropoda", "Velociraptor")
    )
    lins[[taxon]]
  })
  result <- focal_distances("Tyrannosaurus", "Velociraptor", progress = FALSE)
  expect_no_error(print(result))
  expect_invisible(print(result))
})

# ── Network tests (skipped on CRAN) ──────────────────────────────────────────

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

test_that("get_taxonomicon_id skips entry whose lineage has no Biota", {
  clear_cache()
  fake_html <- '
    <html><body><table>
      <tr>
        <td>Carnotaurus - animal - dinosaur</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=12345&src=0">tree</a></td>
      </tr>
    </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)
  mockery::stub(get_taxonomicon_id, "get_lineage_by_id",
                function(...) c("NotBiota", "Animalia", "Dinosauria"))
  result <- get_taxonomicon_id("Carnotaurus")
  expect_null(result)
})

test_that("get_lineage_by_id truncates at own id when present in links", {
  clear_cache()
  fake_html <- '
    <html><body>
      <a href="TaxonTree.aspx?id=1&src=0">Biota</a>
      <a href="TaxonTree.aspx?id=2&src=0">Animalia</a>
      <a href="TaxonTree.aspx?id=99&src=0">Carnotaurus</a>
      <a href="TaxonTree.aspx?id=100&src=0">SomeChild</a>
    </body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) fake_html)
  result <- get_lineage_by_id("99")
  expect_true("Carnotaurus" %in% result)
  expect_false("SomeChild" %in% result)
})

test_that("get_lineage handles binomial taxon name correctly", {
  clear_cache()
  mockery::stub(get_lineage, "get_taxonomicon_id", function(...) "12345")
  mockery::stub(get_lineage, "get_lineage_by_id",
                function(...) c("Biota", "Animalia", "Dinosauria",
                                "Theropoda", "Carnotaurus sastrei"))
  result <- get_lineage("Carnotaurus sastrei")
  expect_true("Carnotaurus sastrei" %in% result)
})

test_that("get_lineage returns NULL when get_lineage_by_id returns NULL", {
  clear_cache()
  mockery::stub(get_lineage, "get_taxonomicon_id", function(...) "12345")
  mockery::stub(get_lineage, "get_lineage_by_id", function(...) NULL)
  result <- get_lineage("Carnotaurus")
  expect_null(result)
})

test_that("get_lineage returns single-node lineage when lineage_by_id returns empty", {
  clear_cache()
  mockery::stub(get_lineage, "get_taxonomicon_id", function(...) "12345")
  mockery::stub(get_lineage, "get_lineage_by_id", function(...) character(0))
  result <- get_lineage("Carnotaurus")
  expect_equal(result, "Carnotaurus")
})

test_that("get_lineage appends taxon name when not found in scraped lineage", {
  clear_cache()
  mockery::stub(get_lineage, "get_taxonomicon_id", function(...) "12345")
  mockery::stub(get_lineage, "get_lineage_by_id",
                function(...) c("Biota", "Animalia", "Dinosauria", "Theropoda"))
  result <- get_lineage("Carnotaurus")
  expect_equal(tail(result, 1), "Carnotaurus")
})

test_that("get_lineage_by_id returns NULL when all links are filtered out", {
  clear_cache()
  fake_html <- '
    <html><body>
      <a href="TaxonTree.aspx?id=1&src=0">Go to</a>
      <a href="TaxonTree.aspx?id=2&src=0">[unranked]</a>
    </body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) fake_html)
  result <- get_lineage_by_id("99999")
  expect_null(result)
})

test_that("get_taxonomicon_id warns on multiple biological entries (coverage)", {
  clear_cache()
  assign("lin_111", c("Biota", "Animalia", "Fake1"), envir = taxodist:::.taxodist_cache)
  assign("lin_222", c("Biota", "Animalia", "Fake2"), envir = taxodist:::.taxodist_cache)

  fake_html <- '
    <html><body><table>
      <tr>
        <td>Nereis - animal - one</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=111&src=0">tree</a></td>
      </tr>
      <tr>
        <td>Nereis - animal - two</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=222&src=0">tree</a></td>
      </tr>
    </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)
  expect_warning(
    get_taxonomicon_id("Nereis"),
    "Multiple valid biological entries"
  )
})

test_that("deduplication preserves order", {
  clear_cache()
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) {
    structure(list(), class = "response")
  })
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) "x")
  mockery::stub(get_lineage_by_id, "rvest::read_html", function(...) {
    xml2::read_html('
        <html><body>
          <a href="TaxonTree.aspx?id=1&src=0">Biota</a>
          <a href="TaxonTree.aspx?id=2&src=0">Animalia</a>
          <a href="TaxonTree.aspx?id=3&src=0">Uropygi</a>
          <a href="TaxonTree.aspx?id=3&src=0">Uropygi</a>
          <a href="TaxonTree.aspx?id=4&src=0">Thelyphonida</a>
        </body></html>')
  })
  result <- get_lineage_by_id("4")
  expect_equal(result, c("Biota", "Animalia", "Uropygi", "Thelyphonida"))
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

test_that("get_lineage accepts numeric IDs directly without searching", {
  clear_cache()
  mockery::stub(get_lineage, "get_lineage_by_id", function(...) c("Biota", "Bacteria"))
  result <- get_lineage("71320")
  expect_equal(result, c("Biota", "Bacteria"))
})

test_that("get_lineage_by_id returns NULL for non-numeric strings", {
  expect_null(get_lineage_by_id("Bacteria"))
  expect_null(get_lineage_by_id("123x"))
  expect_null(get_lineage_by_id("   "))
})

test_that("get_taxonomicon_id follows taxonomic redirects", {
  clear_cache()
  assign("lin_16197", c("Biota", "Animalia", "Uropygi"), envir = taxodist:::.taxodist_cache)
  fake_html <- '
    <html><body><table>
      <tr>
        <td>Thelyphonida see Uropygi</td>
        <td>
          <a class="Invalid" href="TaxonTree.aspx?id=123&src=0">old</a>
          <a class="Valid" href="TaxonTree.aspx?id=16197&src=0">tree</a>
        </td>
      </tr>
    </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)
  result <- get_taxonomicon_id("Thelyphonida")
  expect_equal(result, "16197")
})

test_that("get_taxonomicon_id skips rows with no Valid links", {
  clear_cache()
  assign("lin_222", c("Biota", "Animalia"), envir = taxodist:::.taxodist_cache)

  fake_html <- '
    <html><body><table>
      <tr>
        <td>Invalid taxon</td>
        <td><a class="Invalid" href="TaxonTree.aspx?id=111&src=0">skip me</a></td>
      </tr>
      <tr>
        <td>Good taxon - animal</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=222&src=0">tree</a></td>
      </tr>
    </table></body></html>'

  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)

  result <- get_taxonomicon_id("Good taxon")
  expect_equal(result, "222")
})

test_that("get_taxonomicon_id skips valid links missing numeric IDs", {
  clear_cache()
  assign("lin_333", c("Biota", "Animalia"), envir = taxodist:::.taxodist_cache)

  fake_html <- '
    <html><body><table>
      <tr>
        <td>Missing ID taxon</td>
        <td><a class="Valid" href="TaxonTree.aspx?wrongparam=abc">skip me</a></td>
      </tr>
      <tr>
        <td>Good taxon - animal</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=333&src=0">tree</a></td>
      </tr>
    </table></body></html>'

  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)

  result <- get_taxonomicon_id("Good taxon")
  expect_equal(result, "333")
})

# ── Fallback handling for NAs and NULLs ───────────────────────────────────────

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

test_that("get_lineage_by_id parses #divPageContent path with nav header and ᵀ marker", {
  clear_cache()

  fake_html <- '
    <html><body>
      <div id="ctl00_divSubject"><b>Drosophila</b></div>
      <div id="divPageContent">
  HierarchyNomenclature
  Classification by:
  Systema Naturae 2000 cactophilic {Drosophila} descriptions
  Natura - nature
  actualia - actual entities
  Clade Biota Wagner 2004
  Kingdom Animalia
  Phylum Arthropoda
  Class Insecta
  Order Diptera
  Family Drosophilidae Rondani, 1856
  Genus Drosophila\u1D40 Fall\u00e9n, 1823
  Drosophila melanogaster Meigen, 1830
      </div>
    </body></html>'

  fake_response <- structure(list(status_code = 200L), class = "response")

  mockery::stub(get_lineage_by_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) fake_html)

  lin <- get_lineage_by_id("28940")

  expect_true("Biota" %in% lin)
  expect_true("Drosophila" %in% lin)
  expect_false(any(grepl("melanogaster", lin)))
})
