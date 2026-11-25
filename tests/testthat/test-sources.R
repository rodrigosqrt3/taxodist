# Open Tree Taxonomy and NCBI Taxonomy sources. Requests are answered from
# real responses stored in fixtures/sources (see tools/record_source_fixtures.R);
# a request without a stored response behaves as a failed request.

source_fixture <- function(...) test_path("fixtures", "sources", ...)

read_fixture_response <- function(file) {
  path <- source_fixture(file)
  if (!file.exists(path)) return(list(status = NA_integer_, text = NULL))
  list(
    status = 200L,
    text = paste(readLines(path, encoding = "UTF-8", warn = FALSE), collapse = "\n")
  )
}

fixture_request <- function(source, method, url, body = NULL) {
  if (grepl("/taxonomy/about$", url)) return(read_fixture_response("ott_about.json"))
  if (grepl("/tnrs/match_names$", url)) return(read_fixture_response("ott_match_names.json"))
  if (grepl("/taxonomy/taxon_info$", url)) {
    return(read_fixture_response(paste0("ott_taxon_info_", body$ott_id, ".json")))
  }
  if (grepl("esearch.fcgi", url, fixed = TRUE)) {
    term <- utils::URLdecode(sub(".*[?&]term=([^&]*).*", "\\1", url))
    name <- tolower(sub("^\"(.*)\"\\[.*$", "\\1", term))
    field <- if (grepl("[Scientific Name]", term, fixed = TRUE)) "scientific" else "all"
    return(read_fixture_response(
      paste0("ncbi_esearch_", gsub("[^a-z0-9]+", "_", name), "_", field, ".json")
    ))
  }
  if (grepl("efetch.fcgi", url, fixed = TRUE)) {
    ids <- strsplit(utils::URLdecode(sub(".*[?&]id=([^&]*).*", "\\1", url)), ",")[[1]]
    taxa <- vapply(ids, function(id) {
      path <- source_fixture(paste0("ncbi_taxon_", id, ".xml"))
      if (!file.exists(path)) return("")
      doc <- xml2::read_xml(path)
      as.character(xml2::xml_find_first(doc, "/TaxaSet/Taxon"))
    }, character(1))
    return(list(status = 200L, text = paste0("<TaxaSet>", paste(taxa, collapse = ""), "</TaxaSet>")))
  }
  list(status = NA_integer_, text = NULL)
}

failing_request <- function(source, method, url, body = NULL) {
  list(status = NA_integer_, text = NULL)
}

# Answers source requests from the fixtures for the rest of the calling test,
# starting from an empty cache.
local_sources <- function(request = fixture_request, env = parent.frame()) {
  clear_cache()
  testthat::local_mocked_bindings(.source_request = request, .env = env)
}

expected_resolutions <- function() {
  jsonlite::read_json(source_fixture("expected_resolutions.json"))
}

NAMES <- c(
  "Homo", "Pan", "Mus", "Drosophila", "Nereis", "Quercus",
  "Fakeosaurus", "human", "Homo sapiens"
)
IDS <- list(ott = c("ott:770309", "417957"), ncbi = c("ncbi:9606", "9598"))

test_that("resolutions reproduce the recorded expectations (shared fixture)", {
  local_sources()
  expected <- expected_resolutions()
  for (key in names(expected)) {
    source <- sub("_within_metazoa$", "", key)
    within <- if (grepl("within", key)) "Metazoa" else NULL
    res <- suppressWarnings(taxo_resolve(
      c(NAMES, IDS[[source]]), source = source, progress = FALSE,
      ambiguity = "first", within = within
    ))
    exp <- expected[[key]]
    expect_identical(attr(res, "source"), exp$source, info = key)
    expect_identical(attr(res, "source_url"), exp$source_url, info = key)
    for (i in seq_along(exp$records)) {
      rec <- exp$records[[i]]
      info <- paste(key, rec$input)
      expect_identical(res$input[[i]], rec$input, info = info)
      expect_identical(res$status[[i]], rec$status, info = info)
      expect_identical(res$id[[i]], if (is.null(rec$id)) NA_character_ else rec$id, info = info)
      expect_identical(res$n_candidates[[i]], as.integer(rec$n_candidates), info = info)
      expect_identical(
        res$lineage[[i]],
        if (is.null(rec$lineage)) NULL else unlist(rec$lineage),
        info = info
      )
      expect_identical(
        res$ranks[[i]],
        if (is.null(rec$ranks)) NULL else vapply(rec$ranks, function(r) if (is.null(r)) NA_character_ else r, ""),
        info = info
      )
      expect_identical(
        res$candidates[[i]]$id,
        vapply(rec$candidates, `[[`, "", "id"),
        info = info
      )
    }
  }
})

test_that("external sources resolve the expected taxa", {
  local_sources()
  ott <- taxo_resolve(c("Homo", "Nereis", "Fakeosaurus"), source = "ott", progress = FALSE)
  expect_identical(ott$id, c("ott:770309", "ott:305738", NA))
  expect_identical(ott$status, c("resolved", "resolved", "unresolved"))
  expect_identical(utils::tail(ott$lineage[[1]], 3L), c("Hominidae", "Homininae", "Homo"))
  expect_identical(ott$lineage[[1]][[1]], "life")
  expect_identical(utils::tail(ott$ranks[[1]], 1L), "Genus")
  expect_match(attr(ott, "source"), "^Open Tree Taxonomy 3\\.7")

  ncbi <- suppressWarnings(taxo_resolve(c("Mus", "Drosophila", "human"), source = "ncbi", progress = FALSE))
  # The nominotypical subgenus Mus (Mus) is not counted as a homonym.
  expect_identical(ncbi$id[[1]], "ncbi:10088")
  expect_identical(ncbi$status[[1]], "resolved")
  # Drosophila is a true homonym (a fly and a fungus); the fly has the lower ID.
  expect_identical(ncbi$id[[2]], "ncbi:7215")
  expect_identical(ncbi$status[[2]], "ambiguous")
  # A common name is matched through "All Names".
  expect_identical(ncbi$resolved_name[[3]], "Homo sapiens")
  expect_identical(ncbi$lineage[[1]][[1]], "cellular organisms")
  expect_identical(attr(ncbi, "source"), "NCBI Taxonomy")
})

test_that("resolutions from other sources feed every matrix function", {
  local_sources()
  res <- taxo_resolve(c("Homo", "Pan", "Mus", "Quercus"), source = "ott", progress = FALSE)
  d <- distance_matrix(res, progress = FALSE)
  expect_s3_class(d, "dist")
  m <- as.matrix(d)
  expect_lt(m["Homo", "Pan"], m["Homo", "Mus"])
  expect_lt(m["Homo", "Mus"], m["Homo", "Quercus"])
  expect_s3_class(taxo_tree(res, progress = FALSE), "phylo")

  bundle <- taxo_bundle(res, progress = FALSE)
  expect_match(bundle$source$name, "^Open Tree Taxonomy")
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  write_taxodist_bundle(bundle, file)
  back <- read_taxodist_bundle(file)
  expect_identical(back$source$name, bundle$source$name)
  expect_identical(back$resolution$id, res$id)
  expect_identical(as.vector(back$matrix), as.vector(d))

  bundle2 <- taxo_bundle(c("Homo", "Pan"), source = "ncbi", progress = FALSE)
  expect_identical(bundle2$source$name, "NCBI Taxonomy")
})

test_that("failed requests are retrieval errors", {
  local_sources(failing_request)
  ott <- taxo_resolve(c("Homo", "770309"), source = "ott", progress = FALSE)
  expect_identical(ott$status, c("retrieval_error", "retrieval_error"))
  expect_identical(attr(ott, "source"), "Open Tree Taxonomy")
  ncbi <- taxo_resolve("Homo", source = "ncbi", progress = FALSE)
  expect_identical(ncbi$status, "retrieval_error")
})

test_that("a candidate whose lineage cannot be retrieved is a retrieval error", {
  local_sources(function(source, method, url, body = NULL) {
    if (grepl("taxon_info", url)) return(list(status = 500L, text = ""))
    fixture_request(source, method, url, body)
  })
  res <- taxo_resolve("Homo", source = "ott", progress = FALSE)
  expect_identical(res$status, "retrieval_error")
  expect_identical(res$candidates[[1]]$id, "ott:770309")
  inside <- taxo_resolve("Homo", source = "ott", progress = FALSE, within = "Metazoa")
  expect_identical(inside$status, "retrieval_error")
})

test_that("within restricts candidates of external sources", {
  local_sources()
  plants <- taxo_resolve(c("Homo", "Quercus"), source = "ott", progress = FALSE,
                         within = "Chloroplastida")
  expect_identical(plants$status, c("unresolved", "resolved"))
  expect_identical(plants$n_candidates, c(0L, 1L))
})

test_that("taxo_sources describes every source", {
  s <- taxo_sources()
  expect_identical(s$id, c("taxonomicon", "ott", "ncbi"))
  expect_true(all(nzchar(s$citation)))
  expect_named(s, c("id", "name", "url", "id_prefix", "ranks", "licence", "citation"))
})

test_that("sources are validated", {
  expect_error(taxo_resolve("Homo", source = "gbif"), "Unknown taxonomic source")
  expect_error(taxo_resolve("Homo", source = c("ott", "ncbi")), "single string")
  expect_error(taxo_resolve("Homo", source = NA_character_), "single string")
})

test_that("source ranks follow the Taxonomicon convention", {
  expect_identical(
    taxodist:::.normalize_source_rank(c("genus", "no rank", "clade", "cellular root",
                                        "", NA, "Superfamily", "no rank - terminal")),
    c("Genus", NA, "Clade", NA, NA, NA, "Superfamily", NA)
  )
})

test_that("nested candidates of the same name are not homonyms", {
  outer <- taxodist:::.outermost_candidates
  expect_identical(
    outer(list(c("a", "Mus"), c("a", "Mus", "Mus (Mus)"), c("b", "Mus"), NULL)),
    c(TRUE, FALSE, TRUE, TRUE)
  )
  expect_identical(outer(list(c("a", "x"))), TRUE)
})

test_that("the Taxonomicon ignores nominotypical subgenera when resolving", {
  resolve <- taxo_resolve
  mockery::stub(resolve, ".taxo_search_details", function(...) list(
    status = "ok",
    results = data.frame(id = c("2", "1"), name = c("Mus (Mus)", "Mus"),
                         stringsAsFactors = FALSE)
  ))
  mockery::stub(resolve, "get_lineage_by_id", function(id, ...) {
    if (id == "1") c("Biota", "Rodentia", "Mus") else c("Biota", "Rodentia", "Mus", "Mus (Mus)")
  })
  mockery::stub(resolve, ".cached_ranks", function(id, lineage) rep(NA_character_, length(lineage)))
  res <- resolve("Mus", progress = FALSE)
  expect_identical(res$status, "resolved")
  expect_identical(res$id, "1")
  expect_identical(res$n_candidates, 1L)
})

test_that("the NCBI request rate depends on the API key", {
  old <- options(taxodist.ncbi_api_key = "", taxodist.min_interval = 0)
  on.exit(options(old), add = TRUE)
  expect_equal(taxodist:::.source_interval("ncbi"), 0.4)
  options(taxodist.ncbi_api_key = "abc")
  expect_equal(taxodist:::.source_interval("ncbi"), 0.12)
  expect_match(taxodist:::.ncbi_url("esearch.fcgi", c(db = "taxonomy")), "api_key=abc")
  expect_equal(taxodist:::.source_interval("ott"), 0)
  options(taxodist.min_interval = "bad")
  expect_equal(taxodist:::.source_interval("ott"), 0)
  options(taxodist.ncbi_api_key = NA)
  expect_identical(taxodist:::.ncbi_api_key(), "")
})

test_that("NCBI merged identifiers and malformed responses are handled", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  taxodist:::.ncbi_parse_taxa("not xml")
  expect_null(taxodist:::.source_cache_get("src_ncbi_lin_1"))
  xml <- paste0(
    "<TaxaSet><Taxon><TaxId>2</TaxId><ScientificName>B</ScientificName>",
    "<Rank>genus</Rank><AkaTaxIds><TaxId>3</TaxId></AkaTaxIds>",
    "<LineageEx><Taxon><TaxId>1</TaxId><ScientificName>A</ScientificName>",
    "<Rank>no rank</Rank></Taxon></LineageEx></Taxon></TaxaSet>"
  )
  taxodist:::.ncbi_parse_taxa(xml)
  merged <- taxodist:::.source_cache_get("src_ncbi_lin_3")
  expect_identical(as.character(merged), c("A", "B"))
  expect_identical(attr(merged, "ranks"), c(NA, "Genus"))
  expect_null(taxodist:::.parse_json(NULL))
  expect_null(taxodist:::.parse_json("{not json"))
})

test_that("external requests are built and throttled per source", {
  seen <- character(0)
  captured <- NULL
  request <- taxodist:::.source_request
  mockery::stub(request, ".throttle", function(source, interval) seen <<- c(seen, source))
  mockery::stub(request, "do.call", function(what, args) {
    captured <<- args
    structure(
      list(
        url = "http://example.org", status_code = 200L,
        headers = list(`content-type` = "application/json; charset=utf-8"),
        content = charToRaw("{\"a\":1}")
      ),
      class = "response"
    )
  })
  res <- request("ott", "POST", "http://example.org", body = list(a = 1))
  expect_identical(res, list(status = 200L, text = "{\"a\":1}"))
  expect_identical(seen, "ott")
  expect_identical(as.character(captured$body), "{\"a\":1}")
  expect_identical(captured[[1]], "POST")
  # 429 (too many requests) is retried; other client errors are not.
  expect_false(429L %in% captured$terminate_on)
  expect_true(all(c(400L, 404L) %in% captured$terminate_on))

  mockery::stub(request, "do.call", function(...) stop("offline"))
  expect_identical(request("ncbi", "GET", "http://example.org"),
                   list(status = NA_integer_, text = NULL))
})

test_that("compare_sources compares distances across sources", {
  local_sources()
  set.seed(42)
  cmp <- compare_sources(
    c("Homo", "Pan", "Mus", "Quercus", "Fakeosaurus", "Homo"),
    sources = c("ott", "ncbi"), permutations = 99
  )
  expect_s3_class(cmp, "taxodist_source_comparison")
  expect_identical(cmp$shared, c("Homo", "Pan", "Mus", "Quercus"))
  expect_identical(nrow(cmp$taxa), 6L)
  expect_true(all(c("status_ott", "id_ncbi", "depth_ncbi") %in% names(cmp$taxa)))
  expect_identical(cmp$agreement$n_taxa, 4L)
  expect_true(cmp$agreement$spearman > 0.9)
  expect_true(cmp$agreement$p_value > 0 && cmp$agreement$p_value <= 1)
  expect_output(print(cmp), "ott")
  expect_message(print(cmp), "Taxa resolved in every source: 4 of 6")

  no_test <- compare_sources(c("Homo", "Pan", "Mus"), sources = c("ott", "ncbi"),
                             permutations = 0, within = c(ott = "Metazoa"), progress = FALSE)
  expect_true(is.na(no_test$agreement$p_value))
  shared_within <- compare_sources(c("Homo", "Pan", "Mus"), sources = c("ott", "ncbi"),
                                   permutations = 0, within = "Metazoa", progress = FALSE)
  expect_identical(shared_within$shared, c("Homo", "Pan", "Mus"))
})

test_that("compare_sources validates its arguments", {
  expect_error(compare_sources(c("a", "b", "c"), sources = "ott"), "at least two")
  expect_error(compare_sources(c("a", "b", "c"), sources = c("ott", "ott")), "at least two")
  expect_error(compare_sources(c("a", "b", "c"), sources = c("ott", "gbif")), "Unknown")
  expect_error(compare_sources(c("a", "b"), sources = c("ott", "ncbi")), "at least three")
  expect_error(compare_sources(c("a", "b", "c"), sources = c("ott", "ncbi"), permutations = -1),
               "non-negative")
  local_sources()
  expect_error(
    compare_sources(c("Homo", "Pan", "Mus"), sources = c("ott", "ncbi"),
                    within = c("a", "b"), progress = FALSE),
    "named by source"
  )
})

test_that("the Mantel statistic handles degenerate inputs", {
  mantel <- taxodist:::.mantel
  a <- stats::as.dist(matrix(c(0, 1, 2, 1, 0, 3, 2, 3, 0), 3))
  b <- stats::as.dist(matrix(c(0, Inf, 2, Inf, 0, 3, 2, 3, 0), 3))
  expect_identical(mantel(a, b, 10L), list(r = NA_real_, p = NA_real_))
  flat <- stats::as.dist(matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3))
  expect_true(is.na(mantel(a, flat, 10L)$r))
  set.seed(1)
  expect_equal(mantel(a, a, 9L)$r, 1)
})

test_that("verbose resolution reports each request", {
  local_sources()
  expect_message(
    taxo_resolve(c("Homo", "Pan"), source = "ott", progress = FALSE, verbose = TRUE),
    "Matching 2 names in the Open Tree Taxonomy"
  )
  expect_message(
    taxo_resolve("Mus", source = "ott", progress = FALSE, verbose = TRUE),
    "Retrieving lineage of ott:1068778"
  )
  expect_message(
    taxo_resolve("Quercus", source = "ncbi", progress = FALSE, verbose = TRUE),
    "Searching NCBI Taxonomy for Quercus"
  )
  clear_cache()
  expect_message(
    taxo_resolve("Quercus", source = "ncbi", progress = FALSE, verbose = TRUE),
    "Retrieving 1 lineage from NCBI Taxonomy"
  )
})

test_that("missing metadata and failed lineage batches are tolerated", {
  local_sources(function(source, method, url, body = NULL) {
    if (grepl("/taxonomy/about$", url)) {
      return(list(status = 200L, text = "{\"version\": \"9.9\", \"source\": \"ott9.9\"}"))
    }
    if (grepl("efetch.fcgi", url, fixed = TRUE)) return(list(status = 503L, text = ""))
    fixture_request(source, method, url, body)
  })
  about <- taxodist:::.ott_about()
  expect_identical(about$name, "Open Tree Taxonomy 9.9 (ott9.9)")
  expect_identical(about$url, "https://tree.opentreeoflife.org/about/taxonomy-version")
  # The name was found but its lineage could not be retrieved.
  res <- taxo_resolve("Quercus", source = "ncbi", progress = FALSE)
  expect_identical(res$status, "retrieval_error")
})
