# Real TaxonTree and TaxonList pages, trimmed to the parts the parser reads,
# are stored in tests/testthat/fixtures. expected_lineages.json holds the raw
# lineage (before trimming at Biota) expected from each tree page; the same
# files are used by the Python and Julia implementations to keep the three
# parsers identical.

fixture <- function(name) test_path("fixtures", name)

read_fixture_page <- function(name) {
  xml2::read_html(fixture(name), encoding = "UTF-8")
}

test_that("structured parser reproduces the expected lineage of every fixture", {
  expected <- jsonlite::read_json(fixture("expected_lineages.json"), simplifyVector = TRUE)
  expect_gt(length(expected), 5L)
  for (id in names(expected)) {
    page <- read_fixture_page(paste0("tree_", id, ".html"))
    expect_identical(
      as.character(taxodist:::.lineage_from_page(page, id)),
      expected[[id]],
      info = paste("fixture", id)
    )
  }
})

test_that("structured parser keeps clean names despite authors and markers", {
  passer <- taxodist:::.lineage_from_page(read_fixture_page("tree_59107.html"), "59107")
  expect_true(all(c("Passeres", "Passerida", "Terrestrornithes") %in% passer))
  expect_false(any(grepl("[ ,(]", passer[passer != "Passer"])))
  expect_identical(utils::tail(passer, 1L), "Passer")
  # Children of the page's taxon (species) are excluded.
  expect_false(any(startsWith(passer, "Passer ")))

  microtus <- taxodist:::.lineage_from_page(read_fixture_page("tree_62669.html"), "62669")
  expect_true("Supramyomorpha" %in% microtus)
  expect_identical(utils::tail(microtus, 1L), "Microtus")

  lipotes <- taxodist:::.lineage_from_page(read_fixture_page("tree_68670.html"), "68670")
  expect_identical(utils::tail(lipotes, 3L), c("Lipotoidea", "Lipotidae", "Lipotes"))

  anamorph <- taxodist:::.lineage_from_page(read_fixture_page("tree_204388.html"), "204388")
  expect_identical(utils::tail(anamorph, 1L), "Epidermophyton")
})

test_that("structured parser appends the subject when its own row is missing", {
  yi <- taxodist:::.lineage_from_page(read_fixture_page("tree_7248674.html"), "7248674")
  expect_identical(utils::tail(yi, 2L), c("Scansoriopterygidae", "Yi"))
  expect_false("[plesion]" %in% yi)
})

test_that("structured parser handles species pages", {
  homo <- taxodist:::.lineage_from_page(read_fixture_page("tree_66295.html"), "66295")
  expect_identical(utils::tail(homo, 2L), c("Homo", "Homo sapiens"))
})

test_that("get_lineage_by_id uses the structured parser and caches both variants", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  html <- paste(readLines(fixture("tree_59107.html"), encoding = "UTF-8"), collapse = "\n")
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) structure(list(), class = "response"))
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) html)
  cleaned <- get_lineage_by_id("59107")
  raw <- get_lineage_by_id("59107", clean = FALSE)
  expect_identical(cleaned[[1]], "Biota")
  expect_identical(raw[[1]], "Natura")
  expect_identical(utils::tail(raw, length(cleaned)), cleaned)
})

test_that(".finalize_lineage keeps name-ancestors of species and subspecies", {
  finalize <- taxodist:::.finalize_lineage
  expect_identical(
    finalize(c("Biota", "Homo", "Homo sapiens"), "Homo sapiens sapiens"),
    c("Biota", "Homo", "Homo sapiens", "Homo sapiens sapiens")
  )
  expect_identical(
    finalize(c("Biota", "Microtus", "Microtus (Microtus)", "Microtus arvalis"), "Microtus arvalis"),
    c("Biota", "Microtus", "Microtus (Microtus)", "Microtus arvalis")
  )
  # A subgenus is a child of a queried genus and is dropped.
  expect_identical(
    finalize(c("Biota", "Microtus", "Microtus (Microtus)"), "Microtus"),
    c("Biota", "Microtus")
  )
  # Species of another genus are not name-ancestors.
  expect_identical(
    finalize(c("Biota", "Homo", "Pan troglodytes"), "Homo sapiens"),
    c("Biota", "Homo", "Homo sapiens")
  )
})

test_that("fallback text parser extracts names regardless of author format", {
  extract <- taxodist:::.extract_scientific_name
  raw <- c(
    "Ichthyosauria de Blainville, 1835", "Foraminifera d'Orbigny, 1826",
    "Xenus Özdikmen, 2009", "Xenus O'Brien, 1990",
    "Xidae van der Hoeven, 1850", "Homo sapiens Linnaeus, 1758",
    "Homo sapiens sapiens Linnaeus, 1758",
    "Alismatales R.Br. ex Bercht. & J.Presl", "Firmicutes corrig. Gibbons",
    "Scombriformes - mackerels", "Romeriida Gauthier,",
    "Accipitres (C. Linnaeus, 1758)", "neomura Cavalier-Smith, 1987",
    "Pan-Aves Gauthier & de Queiroz, 2001"
  )
  expect_identical(extract(raw), c(
    "Ichthyosauria", "Foraminifera", "Xenus", "Xenus", "Xidae",
    "Homo sapiens", "Homo sapiens sapiens", "Alismatales", "Firmicutes",
    "Scombriformes", "Romeriida", "Accipitres", "neomura", "Pan-Aves"
  ))
  # Quoted and bracketed lines are left for the caller's filters.
  expect_identical(extract(c('"Setaphyta"', "[unranked]")), c('"Setaphyta"', "[unranked]"))
  expect_identical(extract(NA_character_), NA_character_)
})

test_that("fallback text parser strips double daggers and qualified markers", {
  fake_html <- paste0(
    '<html><body><div id="ctl00_divSubject"><b>Lipotes</b></div>',
    '<div id="divPageContent">Natura - nature\nClade Biota Wagner 2004\n',
    'Superfamily ‡?Lipotoidea de Muizon, 1988\n',
    '[plesion] Family †Lipotidae Zhou, 1978\n',
    'Genus ‡?Lipotesᵀ Miller, 1918\n</div></body></html>'
  )
  page <- xml2::read_html(fake_html)
  expect_identical(
    as.character(taxodist:::.lineage_from_page(page, "1")),
    c("Natura", "Biota", "Lipotoidea", "Lipotidae", "Lipotes")
  )
})

test_that("regexEscape quotes every metacharacter", {
  for (s in c("Aus (Bus) cus", "a.b", "x|y", "a[b]", "a{2}", "a+b", "a^b$", "a\\b", "a*?")) {
    expect_true(grepl(paste0("^", taxodist:::regexEscape(s), "$"), s), info = s)
  }
  expect_true(taxodist:::.contains_word("Microtus (Agricola)", "Microtus (Agricola)"))
  expect_true(taxodist:::.contains_word("Aus (Bus) cus", "Aus (Bus)"))
  expect_false(taxodist:::.contains_word("Homopus", "Homo"))
  expect_true(taxodist:::.contains_word("Émys", "Émys"))
  expect_false(taxodist:::.contains_word("xémys", "mys"))
})

test_that("taxo_resolve reports non-biological candidates as unresolved", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  assign("lin_1", c("Natura", "Mundus", "Venus"), envir = taxodist:::.taxodist_cache)
  mockery::stub(taxo_resolve, ".taxo_search_details", function(...) {
    list(status = "ok", results = data.frame(id = "1", name = "Venus - planet"))
  })
  result <- taxo_resolve("Venus", progress = FALSE)
  expect_identical(result$status, "unresolved")
  expect_identical(result$n_candidates, 0L)

  # The record must survive bundle validation and a JSON round trip.
  ok <- taxo_from_lineages(list(Homo = c("Biota", "Homo"), Pan = c("Biota", "Pan")))
  mixed <- rbind(ok, result)
  attr(mixed, "source") <- "test"
  bundle <- taxo_bundle(mixed, progress = FALSE)
  tmp <- tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  write_taxodist_bundle(bundle, tmp)
  expect_identical(read_taxodist_bundle(tmp)$resolution$status[[3]], "unresolved")
})

test_that("taxo_resolve keeps retrieval_error when a candidate cannot be fetched", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  mockery::stub(taxo_resolve, ".taxo_search_details", function(...) {
    list(status = "ok", results = data.frame(id = "abc", name = "Broken"))
  })
  result <- taxo_resolve("Broken", progress = FALSE)
  expect_identical(result$status, "retrieval_error")
})

test_that("invalid or rejected homonyms are ranked after valid ones", {
  order_fn <- taxodist:::.candidate_order
  desc <- c(
    "Amanita (Genus) Dill., nom. inval. see Amanita Pers.",
    "Amanita (Genus) Pers. (1797), nom. cons. (fungus)",
    "Amanita (Genus) Foo, nom. rej."
  )
  expect_identical(order_fn(desc), c(2L, 1L, 3L))
  expect_identical(order_fn(c("a", "b")), c(1L, 2L))
})

test_that("cache files from an older parser keep only taxon IDs", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(list(id_Passer = "59107", lin_59107 = c("Biota", "Passeres (C.")), tmp)
  expect_warning(load_cache(tmp), "older version")
  expect_true(exists("id_Passer", envir = taxodist:::.taxodist_cache))
  expect_false(exists("lin_59107", envir = taxodist:::.taxodist_cache))

  clear_cache()
  saveRDS(list(), tmp)
  expect_no_warning(load_cache(tmp))

  assign("lin_1", c("Biota", "A"), envir = taxodist:::.taxodist_cache)
  save_cache(tmp)
  clear_cache()
  expect_no_warning(load_cache(tmp))
  expect_identical(get("lin_1", envir = taxodist:::.taxodist_cache), c("Biota", "A"))
  expect_false(exists("meta_parser_version", envir = taxodist:::.taxodist_cache))
})

test_that("requests are throttled by taxodist.min_interval", {
  throttle <- taxodist:::.throttle
  old <- options(taxodist.min_interval = 0.2)
  on.exit(options(old), add = TRUE)
  throttle()
  start <- Sys.time()
  throttle()
  expect_gte(as.numeric(difftime(Sys.time(), start, units = "secs")), 0.15)

  options(taxodist.min_interval = 0)
  start <- Sys.time()
  throttle()
  expect_lt(as.numeric(difftime(Sys.time(), start, units = "secs")), 0.1)
})
