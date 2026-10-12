# Retrieval and parsing of The Taxonomicon pages (mocked)

test_that("get_taxonomicon_id returns NULL on network failure", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
  clear_cache()
  mockery::stub(lookup_id, "httr::RETRY", function(...) NULL)
  expect_warning(result <- get_id("Tyrannosaurus"), "Cannot reach")
  expect_null(result)
})

test_that("get_taxonomicon_id returns NULL on bad status", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 404L)
  expect_warning(result <- get_id("Tyrannosaurus"), "Cannot reach")
  expect_null(result)
})

test_that("get_taxonomicon_id returns NULL on an unparseable response", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) "not html")
  mockery::stub(
    get_taxonomicon_id,
    "rvest::read_html",
    function(...) stop("parse failure")
  )
  expect_warning(
    result <- get_id("Tyrannosaurus"),
    "Could not parse"
  )
  expect_null(result)
})

test_that("get_lineage_by_id returns NULL on network failure", {
  clear_cache()
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) NULL)
  result <- get_lineage_by_id("12345")
  expect_null(result)
})

test_that("get_lineage returns NULL when id not found", {
  clear_cache()
  mockery::stub(get_lineage, "get_taxonomicon_id", function(...) NULL)
  result <- get_lineage("Fakeosaurus")
  expect_null(result)
})

test_that("get_taxonomicon_id verbose prints warning on network failure", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
  clear_cache()
  mockery::stub(lookup_id, "httr::RETRY", function(...) NULL)
  expect_warning(get_id("Drosophila", verbose = TRUE), "Cannot reach")
})

test_that("get_taxonomicon_id verbose prints warning on bad status", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 503L)
  expect_warning(get_id("Drosophila", verbose = TRUE), "Cannot reach")
})

test_that("get_lineage_by_id verbose prints warning on network failure", {
  clear_cache()
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) NULL)
  expect_no_error(get_lineage_by_id("00000", verbose = TRUE))
})

test_that("get_lineage_by_id returns NULL when lineage is empty after cleaning", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) fake_response)
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
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) fake_response)
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
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 404L)
  result <- get_lineage_by_id("99999")
  expect_null(result)
})

test_that("get_lineage_by_id verbose prints warning on bad status", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 503L)
  expect_no_error(get_lineage_by_id("99999", verbose = TRUE))
})

test_that("get_lineage_by_id with clean = FALSE returns NULL on network failure", {
  clear_cache()
  assign("lin_clean_test", NULL, envir = taxodist:::.taxodist_cache)
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) NULL)
  result <- get_lineage_by_id("clean_test", clean = FALSE)
  expect_null(result)
})

test_that("get_taxonomicon_id parses HTML and returns id", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
  clear_cache()
  fake_html <- '
    <html><body><table>
      <tr>
        <td>Carnotaurus - animal - dinosaur</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=12345&src=0">tree</a></td>
      </tr>
    </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)
  mockery::stub(lookup_id, "get_lineage_by_id", function(...) c("Biota", "Animalia"))
  result <- get_id("Carnotaurus")
  expect_equal(result, "12345")
})

test_that("get_taxonomicon_id skips astronomical entries", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
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
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)
  mockery::stub(lookup_id, "get_lineage_by_id", function(...) c("Biota", "Animalia"))
  result <- get_id("Carnotaurus", verbose = TRUE)
  expect_equal(result, "12345")
})

test_that("get_taxonomicon_id skips row matching both bio and astronomical keywords", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
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
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)
  mockery::stub(lookup_id, "get_lineage_by_id", function(...) c("Biota", "Animalia"))
  result <- get_id("Pterodactylus", verbose = TRUE)
  expect_equal(result, "42042")
})

test_that("get_taxonomicon_id returns NULL when bio row has no tree link", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
  clear_cache()
  fake_html <- '
    <html><body><table>
      <tr>
        <td>Quercus - plant</td>
        <td><a href="SomeOtherPage.aspx?id=999">no tree link</a></td>
      </tr>
    </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)
  result <- get_id("Quercus", verbose = TRUE)
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
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) fake_response)
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
  mockery::stub(search_details, "httr::RETRY", function(...) stop("Network error"))
  expect_equal(
    search_details("Bacteria", verbose = TRUE)$status,
    "retrieval_error"
  )

  fake_response <- structure(list(), class = "response")
  mockery::stub(search_details, "httr::RETRY", function(...) fake_response)
  mockery::stub(search_details, "httr::status_code", function(...) 503L)
  expect_equal(search_details("Bacteria")$status, "retrieval_error")

  mockery::stub(search_details, "httr::status_code", function(...) 200L)
  mockery::stub(search_details, "httr::content", function(...) "not html")
  mockery::stub(search_details, "rvest::read_html", function(...) {
    stop("parse failure")
  })
  expect_message(
    details <- search_details("Bacteria", verbose = TRUE),
    "Could not parse"
  )
  expect_equal(details$status, "retrieval_error")
})

test_that("taxo_search returns NULL when no matches are found", {
  clear_cache()
  fake_response <- structure(list(), class = "response")
  search_details <- taxodist:::.taxo_search_details
  mockery::stub(search_details, "httr::RETRY", function(...) fake_response)
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
  mockery::stub(search_details, "httr::RETRY", function(...) fake_response)
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

test_that("get_taxonomicon_id skips entry whose lineage has no Biota", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
  clear_cache()
  fake_html <- '
    <html><body><table>
      <tr>
        <td>Carnotaurus - animal - dinosaur</td>
        <td><a class="Valid" href="TaxonTree.aspx?id=12345&src=0">tree</a></td>
      </tr>
    </table></body></html>'
  fake_response <- structure(list(), class = "response")
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)
  mockery::stub(lookup_id, "get_lineage_by_id",
                function(...) c("NotBiota", "Animalia", "Dinosauria"))
  result <- get_id("Carnotaurus")
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
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) fake_response)
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
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) fake_html)
  result <- get_lineage_by_id("99999")
  expect_null(result)
})

test_that("get_taxonomicon_id warns on multiple biological entries (coverage)", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
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
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)
  expect_warning(
    get_id("Nereis"),
    "Multiple valid biological entries"
  )
})

test_that("deduplication preserves order", {
  clear_cache()
  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) {
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
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
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
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)
  result <- get_id("Thelyphonida")
  expect_equal(result, "16197")
})

test_that("get_taxonomicon_id skips rows with no Valid links", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
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
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)

  result <- get_id("Good taxon")
  expect_equal(result, "222")
})

test_that("get_taxonomicon_id skips valid links missing numeric IDs", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
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
  mockery::stub(lookup_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(lookup_id, "httr::status_code", function(...) 200L)
  mockery::stub(lookup_id, "httr::content", function(...) fake_html)

  result <- get_id("Good taxon")
  expect_equal(result, "333")
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

  mockery::stub(get_lineage_by_id, "httr::RETRY", function(...) fake_response)
  mockery::stub(get_lineage_by_id, "httr::status_code", function(...) 200L)
  mockery::stub(get_lineage_by_id, "httr::content", function(...) fake_html)

  lin <- get_lineage_by_id("28940")

  expect_true("Biota" %in% lin)
  expect_true("Drosophila" %in% lin)
  expect_false(any(grepl("melanogaster", lin)))
})

