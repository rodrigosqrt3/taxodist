# Distance definition and distance-based queries

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
      taxodist:::.pairwise_distances(case),
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
  mockery::stub(check_coverage, ".lookup_taxonomicon_id", function(taxon, ...) {
    switch(taxon,
      Fakeosaurus = list(id = NULL, status = "not_found"),
      Offline = list(id = NULL, status = "retrieval_error"),
      list(id = "12345", status = "found")
    )
  })
  result <- check_coverage(c("Tyrannosaurus", "Fakeosaurus", "Offline"))
  expect_type(result, "logical")
  expect_true(result[["Tyrannosaurus"]])
  expect_false(result[["Fakeosaurus"]])
  expect_true(is.na(result[["Offline"]]))
  expect_named(result, c("Tyrannosaurus", "Fakeosaurus", "Offline"))
  expect_equal(check_coverage(character(0)), stats::setNames(logical(0), character(0)))
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

test_that("distance_matrix rejects a malformed taxodist_resolution", {
  invalid <- structure(
    data.frame(input = "Alpha", status = "resolved"),
    class = c("taxodist_resolution", "data.frame")
  )
  expect_error(
    distance_matrix(invalid, progress = FALSE),
    "Invalid.*taxodist_resolution"
  )
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

