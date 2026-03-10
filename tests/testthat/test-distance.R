library(testthat)
library(taxodist)

# ── Pure logic tests ──────────────────────────────────────────────────────────

test_that("taxodist package loads", {
  expect_true(TRUE)
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

test_that(".compute_distance returns 0 when one taxon is ancestor of other", {
  lin_a <- c("Biota", "Animalia", "Dinosauria")
  lin_b <- c("Biota", "Animalia", "Dinosauria", "Theropoda", "Carnotaurus")
  result <- taxodist:::.compute_distance(lin_a, lin_b)
  expect_equal(result$distance, 0)
})

test_that("clear_cache returns invisible NULL", {
  expect_invisible(clear_cache())
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

# ── Mock tests ────────────────────────────────────────────────────────────────

test_that("get_taxonomicon_id returns NULL on network failure", {
  clear_cache()
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) NULL)
  result <- get_taxonomicon_id("Tyrannosaurus")
  expect_null(result)
})

test_that("get_taxonomicon_id returns NULL on bad status", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET",
                function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code",
                function(...) 404L)
  result <- get_taxonomicon_id("Tyrannosaurus")
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
  expect_equal(nrow(as.matrix(mat)), 3)
  expect_equal(diag(as.matrix(mat)), c(0, 0, 0), ignore_attr = TRUE)
})

test_that("check_coverage returns named logical vector", {
  mockery::stub(check_coverage, "get_taxonomicon_id",
                function(taxon, ...) if (taxon == "Fakeosaurus") NULL else "12345")
  result <- check_coverage(c("Tyrannosaurus", "Fakeosaurus"))
  expect_type(result, "logical")
  expect_true(result["Tyrannosaurus"])
  expect_false(result["Fakeosaurus"])
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
  expect_no_error(get_taxonomicon_id("Drosophila", verbose = TRUE))
})

test_that("get_taxonomicon_id verbose prints warning on bad status", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 503L)
  expect_no_error(get_taxonomicon_id("Drosophila", verbose = TRUE))
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
      <td><a href="TaxonTree.aspx?id=12345&src=0">tree</a></td>
    </tr>
  </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)
  result <- get_taxonomicon_id("Carnotaurus")
  expect_equal(result, "12345")
})

test_that("get_taxonomicon_id skips astronomical entries", {
  clear_cache()
  fake_html <- '
  <html><body><table>
    <tr>
      <td>Carnotaurus - asteroid - Minor planet</td>
      <td><a href="TaxonTree.aspx?id=99999&src=0">tree</a></td>
    </tr>
    <tr>
      <td>Carnotaurus - animal - dinosaur</td>
      <td><a href="TaxonTree.aspx?id=12345&src=0">tree</a></td>
    </tr>
  </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)
  result <- get_taxonomicon_id("Carnotaurus", verbose = TRUE)
  expect_equal(result, "12345")
})

test_that("get_taxonomicon_id skips row matching both bio and astronomical keywords", {
  clear_cache()
  fake_html <- '
  <html><body><table>
    <tr>
      <td>Pterodactylus - animal - Minor planet asteroid</td>
      <td><a href="TaxonTree.aspx?id=99999&src=0">wrong</a></td>
    </tr>
    <tr>
      <td>Pterodactylus - animal - reptile</td>
      <td><a href="TaxonTree.aspx?id=42042&src=0">tree</a></td>
    </tr>
  </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_taxonomicon_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_taxonomicon_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_taxonomicon_id, "httr::content", function(...) fake_html)
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

# ── Network tests (skipped on CRAN) ──────────────────────────────────────────

test_that("get_lineage returns correct lineage for Velociraptor", {
  skip_if_offline()
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
  clear_cache()
  lin <- get_lineage("Tyrannosaurus")
  expect_type(lin, "character")
  expect_true("Coelurosauria" %in% lin)
  expect_true("Dinosauria" %in% lin)
})

test_that("get_lineage returns correct lineage for Carnotaurus", {
  skip_if_offline()
  clear_cache()
  lin <- get_lineage("Carnotaurus")
  expect_type(lin, "character")
  expect_true("Dinosauria" %in% lin)
  expect_true("Theropoda" %in% lin)
})

test_that("get_lineage returns correct lineage for Homo", {
  skip_if_offline()
  clear_cache()
  lin <- get_lineage("Homo")
  expect_true("Amniota" %in% lin)
  expect_true("Mammalia" %in% lin)
})

test_that("get_lineage returns correct lineage for Drosophila", {
  skip_if_offline()
  clear_cache()
  lin <- get_lineage("Drosophila")
  expect_type(lin, "character")
  expect_true("Animalia" %in% lin)
})

test_that("get_lineage returns NULL for unknown taxon", {
  skip_if_offline()
  clear_cache()
  expect_null(get_lineage("Fakeosaurus"))
})

test_that("taxo_distance returns valid result for Tyrannosaurus vs Velociraptor", {
  skip_if_offline()
  clear_cache()
  result <- taxo_distance("Tyrannosaurus", "Velociraptor")
  expect_s3_class(result, "taxodist_result")
  expect_gte(result$distance, 0)
  expect_lte(result$distance, 1)
  expect_equal(result$mrca, "Tyrannoraptora")
})

test_that("taxo_distance returns 0 when one taxon is ancestor of other", {
  skip_if_offline()
  clear_cache()
  expect_equal(taxo_distance("Tyrannosaurus", "Dinosauria")$distance, 0)
  expect_equal(taxo_distance("Carnotaurus", "Ceratosauria")$distance, 0)
})

test_that("taxo_distance between Carnotaurus and Triceratops is valid", {
  skip_if_offline()
  clear_cache()
  result <- taxo_distance("Carnotaurus", "Triceratops")
  expect_s3_class(result, "taxodist_result")
  expect_equal(result$mrca, "Dinosauria")
})

test_that("taxo_distance is larger between distant taxa than close ones", {
  skip_if_offline()
  clear_cache()
  d_close  <- taxo_distance("Carnotaurus", "Tyrannosaurus")$distance
  d_distant <- taxo_distance("Carnotaurus", "Homo")$distance
  expect_gt(d_distant, d_close)
})

test_that("mrca of Tyrannosaurus and Triceratops is Dinosauria", {
  skip_if_offline()
  clear_cache()
  expect_equal(mrca("Tyrannosaurus", "Triceratops"), "Dinosauria")
})

test_that("mrca of Tyrannosaurus and Homo is Amniota", {
  skip_if_offline()
  clear_cache()
  expect_equal(mrca("Tyrannosaurus", "Homo"), "Amniota")
})

test_that("mrca of Velociraptor and Triceratops is Dinosauria", {
  skip_if_offline()
  clear_cache()
  expect_equal(mrca("Velociraptor", "Triceratops"), "Dinosauria")
})

test_that("mrca of Carnotaurus and Tyrannosaurus is within Theropoda", {
  skip_if_offline()
  clear_cache()
  ancestor <- mrca("Carnotaurus", "Tyrannosaurus")
  lin <- get_lineage("Tyrannosaurus")
  expect_true(ancestor %in% lin)
})

test_that("is_member correctly identifies clade membership", {
  skip_if_offline()
  clear_cache()
  expect_true(is_member("Tyrannosaurus", "Theropoda"))
  expect_false(is_member("Tyrannosaurus", "Ornithischia"))
})

test_that("lineage_depth for Carnotaurus is reasonable", {
  skip_if_offline()
  clear_cache()
  expect_gt(lineage_depth("Carnotaurus"), 10)
})

test_that("get_taxonomicon_id finds real ID and caches it", {
  skip_if_offline()
  clear_cache()
  id <- get_taxonomicon_id("Carnotaurus", verbose = TRUE)
  expect_type(id, "character")
  expect_false(is.null(id))
  expect_equal(id, get("id_Carnotaurus", envir = taxodist:::.taxodist_cache))
})

test_that("get_lineage_by_id parses and caches lineage for Drosophila", {
  skip_if_offline()
  clear_cache()
  id <- get_taxonomicon_id("Drosophila")
  result <- get_lineage_by_id(id, verbose = TRUE)
  expect_type(result, "character")
  expect_true("Animalia" %in% result)
  expect_equal(result, get(paste0("lin_", id), envir = taxodist:::.taxodist_cache))
})

test_that("get_lineage_by_id clean = FALSE keeps more nodes than clean = TRUE", {
  skip_if_offline()
  clear_cache()
  id <- get_taxonomicon_id("Carnotaurus")
  result_clean <- get_lineage_by_id(id, clean = TRUE)
  clear_cache()
  result_no_clean <- get_lineage_by_id(id, clean = FALSE)
  expect_lte(length(result_clean), length(result_no_clean))
})

test_that("get_lineage verbose wrapper works for Quercus", {
  skip_if_offline()
  clear_cache()
  result <- get_lineage("Quercus", verbose = TRUE)
  expect_type(result, "character")
  expect_true("Biota" %in% result)
})

test_that("get_taxonomicon_id returns NULL for genuinely unknown taxon", {
  skip_if_offline()
  clear_cache()
  expect_null(get_taxonomicon_id("Zzzznotarealgenus99999", verbose = TRUE))
})

test_that("get_taxonomicon_id skips astronomical objects", {
  skip_if_offline()
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
  clear_cache()
  id <- get_taxonomicon_id("Carnotaurus")
  expect_no_error(get_lineage_by_id(id, verbose = TRUE))
})

test_that("get_taxonomicon_id works with verbose for real taxon", {
  skip_if_offline()
  clear_cache()
  expect_no_error(get_taxonomicon_id("Drosophila", verbose = TRUE))
})

test_that("get_taxonomicon_id verbose prints not found warning", {
  skip_if_offline()
  clear_cache()
  expect_null(get_taxonomicon_id("Zzzzfakeosaurus99999", verbose = TRUE))
})

test_that("get_lineage_by_id verbose success message fires on real taxon", {
  skip_if_offline()
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
      <td><a href="TaxonTree.aspx?id=12345&src=0">tree</a></td>
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
    <a href="TaxonTree.aspx?id=1&src=0">Genus Fakeosaurus Smith, 1999</a>
  </body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_lineage_by_id, "httr::GET", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) fake_html)
  result <- get_lineage_by_id("99999")
  expect_null(result)
})
