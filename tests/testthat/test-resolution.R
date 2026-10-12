# ── Name normalisation, lineage truncation, and ID lookup ────────────────────

test_that("internal helpers normalise names and escape cli braces", {
  expect_equal(taxodist:::.normalize_taxon("  Homo   sapiens "), "Homo sapiens")
  expect_equal(taxodist:::.cli_escape("a {b} c"), "a {{b}} c")
  expect_match(
    taxodist:::.search_url("A&B C"),
    "search=A%26B%20C$"
  )
})

test_that("the retry count option is validated", {
  old <- options(taxodist.max_tries = 5)
  on.exit(options(old))
  expect_identical(taxodist:::.max_tries(), 5L)
  options(taxodist.max_tries = "many")
  expect_identical(taxodist:::.max_tries(), 1L)
  options(taxodist.max_tries = 0)
  expect_identical(taxodist:::.max_tries(), 1L)
})

test_that(".finalize_lineage truncates genus queries exactly and case-insensitively", {
  finalize <- taxodist:::.finalize_lineage
  lineage <- c("Biota", "Dinosauria", "[unranked] x", "Tyrannosaurus",
               "Tyrannosaurus rex")
  expect_equal(
    finalize(lineage, "Tyrannosaurus"),
    c("Biota", "Dinosauria", "Tyrannosaurus")
  )
  expect_equal(
    finalize(lineage, "tyrannosaurus"),
    c("Biota", "Dinosauria", "Tyrannosaurus")
  )
  expect_equal(
    finalize(c("Biota", "Blattidae"), "Periplaneta"),
    c("Biota", "Blattidae", "Periplaneta")
  )
})

test_that(".finalize_lineage handles binomial queries", {
  finalize <- taxodist:::.finalize_lineage
  lineage <- c("Biota", "Hominidae", "Homo", "Homo sapiens", "Homo sapiens sapiens")
  expect_equal(
    finalize(lineage, "homo sapiens"),
    c("Biota", "Hominidae", "Homo", "Homo sapiens")
  )
  expect_equal(
    finalize(c("Biota", "Hominidae", "Homo"), "Homo sapiens"),
    c("Biota", "Hominidae", "Homo", "Homo sapiens")
  )
})

test_that("get_lineage normalises case and whitespace before resolving", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  assign("id_Tyrannosaurus", "1", envir = taxodist:::.taxodist_cache)
  assign("id_tyrannosaurus", "1", envir = taxodist:::.taxodist_cache)
  assign("lin_1", c("Biota", "Dinosauria", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)

  expected <- c("Biota", "Dinosauria", "Tyrannosaurus")
  expect_equal(get_lineage("Tyrannosaurus"), expected)
  expect_equal(get_lineage("  Tyrannosaurus  "), expected)
  expect_equal(get_lineage("tyrannosaurus"), expected)
  expect_equal(taxo_distance("Tyrannosaurus", "tyrannosaurus")$distance, 0)
})

test_that("get_lineage validates its input", {
  expect_error(get_lineage(c("A", "B")), "single taxon")
  expect_null(get_lineage(NA_character_))
  expect_null(get_lineage("   "))
})

test_that("get_lineage accepts numeric IDs given as numbers", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  assign("lin_50841", c("Biota", "Tyrannosaurus"), envir = taxodist:::.taxodist_cache)
  expect_equal(get_lineage(50841), c("Biota", "Tyrannosaurus"))
  expect_equal(get_lineage_by_id(50841), c("Biota", "Tyrannosaurus"))
  expect_equal(get_lineage_by_id(" 50841 "), c("Biota", "Tyrannosaurus"))
})

test_that("get_lineage_by_id validates the ID argument", {
  expect_null(get_lineage_by_id(NULL))
  expect_null(get_lineage_by_id(character(0)))
  expect_null(get_lineage_by_id(NA))
  expect_error(get_lineage_by_id(c("1", "2")), "single ID")
})

test_that("cleaned and uncleaned lineages never overwrite each other in the cache", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  fake_html <- paste0(
    '<html><body><div id="ctl00_divSubject"><b>Foo</b></div>',
    '<div id="divPageContent">Natura\nactualia\nBiota\nAnimalia\nFoo\n</div>',
    "</body></html>"
  )
  calls <- 0L
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) {
    calls <<- calls + 1L
    structure(list(), class = "response")
  })
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) fake_html)

  raw <- get_lineage_by_id("999", clean = FALSE)
  cleaned <- get_lineage_by_id("999", clean = TRUE)
  expect_equal(raw, c("Natura", "actualia", "Biota", "Animalia", "Foo"))
  expect_equal(cleaned, c("Biota", "Animalia", "Foo"))
  expect_equal(calls, 1L)

  clear_cache()
  cleaned_first <- get_lineage_by_id("999", clean = TRUE)
  raw_second <- get_lineage_by_id("999", clean = FALSE)
  expect_equal(cleaned_first, cleaned)
  expect_equal(raw_second, raw)
})

test_that("ID lookup reports statuses for every outcome", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  clear_cache()
  on.exit(clear_cache(), add = TRUE)

  mockery::stub(lookup_id, "httr::RETRY", function(...) stop("offline"))
  expect_warning(result <- lookup_id("Alpha"), "Cannot reach")
  expect_equal(result$status, "retrieval_error")
  expect_null(result$id)
})

test_that("ID lookup reports a parse failure as a retrieval error", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  clear_cache()
  mockery::stub(lookup_id, "httr::RETRY", function(...) structure(list(), class = "response"))
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) "x")
  mockery::stub(lookup_id, "rvest::read_html", function(...) stop("bad"))
  expect_warning(result <- lookup_id("Alpha"), "Could not parse")
  expect_equal(result$status, "retrieval_error")
})

test_that("ID lookup separates absent names from failed candidate retrieval", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  fake_html <- '
    <html><body><table>
      <tr><td>Alpha - animal</td>
      <td><a class="Valid" href="TaxonTree.aspx?id=101&src=0">tree</a></td></tr>
    </table></body></html>'
  mockery::stub(lookup_id, "httr::RETRY", function(...) structure(list(), class = "response"))
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)

  mockery::stub(lookup_id, "get_lineage_by_id", function(...) NULL)
  expect_equal(lookup_id("Alpha")$status, "retrieval_error")

  mockery::stub(lookup_id, "get_lineage_by_id", function(...) c("Natura", "Mineralia"))
  expect_equal(lookup_id("Alpha", verbose = TRUE)$status, "not_found")

  mockery::stub(lookup_id, "get_lineage_by_id", function(...) c("Biota", "Alpha"))
  found <- lookup_id("  Alpha ")
  expect_equal(found, list(id = "101", status = "found"))
  expect_equal(get("id_Alpha", envir = taxodist:::.taxodist_cache), "101")
  expect_equal(lookup_id("Alpha")$status, "found")
})

test_that("ID lookup prefers candidates containing the name and escapes braces", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  fake_html <- '
    <html><body><table>
      <tr><td>Nereis {worm} one</td>
      <td><a class="Valid" href="TaxonTree.aspx?id=1&src=0">tree</a></td></tr>
      <tr><td>Nereis {worm} two</td>
      <td><a class="Valid" href="TaxonTree.aspx?id=2&src=0">tree</a></td></tr>
      <tr><td>Other thing</td>
      <td><a class="Valid" href="TaxonTree.aspx?id=3&src=0">tree</a></td></tr>
    </table></body></html>'
  lineages <- list(
    `1` = c("Biota", "Animalia", "Nereis"),
    `2` = c("Biota", "Plantae", "Nereis"),
    `3` = c("Biota", "Fungi", "Other")
  )
  mockery::stub(lookup_id, "httr::RETRY", function(...) structure(list(), class = "response"))
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)
  mockery::stub(lookup_id, "get_lineage_by_id", function(id, ...) lineages[[id]])

  expect_warning(
    result <- lookup_id("Nereis"),
    "Nereis \\{worm\\} one"
  )
  expect_equal(result$id, "1")
})

test_that("taxo_resolve truncates lineages exactly like get_lineage", {
  mockery::stub(taxo_resolve, ".taxo_search_details", function(...) {
    list(status = "ok", results = data.frame(id = "2", name = "Genus Periplaneta"))
  })
  mockery::stub(taxo_resolve, "get_lineage_by_id", function(...) {
    c("Biota", "Insecta", "Blattidae")
  })
  result <- taxo_resolve("Periplaneta", progress = FALSE)
  expect_equal(result$lineage[[1]], c("Biota", "Insecta", "Blattidae", "Periplaneta"))
  expect_equal(result$resolved_name, "Periplaneta")
  expect_equal(result$lineage_depth, 4L)
})

test_that("taxo_resolve resolves normalised duplicates once and keeps raw inputs", {
  searches <- 0L
  mockery::stub(taxo_resolve, ".taxo_search_details", function(taxon, ...) {
    searches <<- searches + 1L
    list(status = "ok", results = data.frame(id = "101", name = taxon))
  })
  mockery::stub(taxo_resolve, "get_lineage_by_id", function(...) c("Biota", "Alpha"))
  result <- taxo_resolve(c("Alpha", " Alpha  "), progress = FALSE)
  expect_equal(searches, 1L)
  expect_equal(result$input, c("Alpha", " Alpha  "))
  expect_equal(result$id, c("101", "101"))
})

test_that("taxo_resolve restricts candidates to a clade with `within`", {
  candidates <- data.frame(
    id = c("201", "202"),
    name = c("Nereis worm", "Nereis butterfly"),
    stringsAsFactors = FALSE
  )
  lineages <- list(
    `201` = c("Biota", "Animalia", "Annelida", "Nereis"),
    `202` = c("Biota", "Plantae", "Nereis")
  )
  mockery::stub(taxo_resolve, ".taxo_search_details", function(...) {
    list(status = "ok", results = candidates)
  })
  mockery::stub(taxo_resolve, "get_lineage_by_id", function(id, ...) lineages[[id]])

  expect_no_warning(
    animal <- taxo_resolve("Nereis", within = "annelida", progress = FALSE)
  )
  expect_equal(animal$status, "resolved")
  expect_equal(animal$id, "201")
  expect_equal(animal$n_candidates, 1L)

  fungus <- taxo_resolve("Nereis", within = "Fungi", progress = FALSE)
  expect_equal(fungus$status, "unresolved")
  expect_equal(fungus$n_candidates, 0L)
  expect_null(fungus$lineage[[1]])

  expect_error(taxo_resolve("Nereis", within = c("A", "B")), "within")
  expect_error(taxo_resolve("Nereis", within = ""), "within")
})

test_that("distance_matrix agrees for names and their taxo_resolve resolution", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  lineages <- list(
    `1` = c("Biota", "Insecta", "Blattidae"),
    `2` = c("Biota", "Insecta", "Blattidae", "Blatta")
  )
  ids <- c(Periplaneta = "1", Blatta = "2")
  for (name in names(ids)) {
    assign(paste0("id_", name), ids[[name]], envir = taxodist:::.taxodist_cache)
    assign(paste0("lin_", ids[[name]]), lineages[[ids[[name]]]],
           envir = taxodist:::.taxodist_cache)
  }
  mockery::stub(taxo_resolve, ".taxo_search_details", function(taxon, ...) {
    list(status = "ok", results = data.frame(id = ids[[taxon]], name = taxon))
  })

  by_name <- distance_matrix(names(ids), progress = FALSE)
  by_resolution <- distance_matrix(
    taxo_resolve(names(ids), progress = FALSE),
    progress = FALSE
  )
  expect_equal(as.matrix(by_name), as.matrix(by_resolution))
  expect_equal(as.matrix(by_name)[1, 2], 1 / 3)
})
