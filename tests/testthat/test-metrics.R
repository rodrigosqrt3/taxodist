# Ranks, alternative distance metrics, schema 1.1 bundles, and taxonomic
# distinctness.

fixture <- function(name) test_path("fixtures", name)

ranked_example <- function() {
  lineages <- list(
    A = c("Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "Felidae", "Panthera", "A"),
    B = c("Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "Felidae", "Panthera", "B"),
    C = c("Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "Canidae", "Canis", "C"),
    D = c("Biota", "Animalia", "Chordata", "Mammalia", "Primates", "Hominidae", "Homo", "D"),
    E = c("Biota", "Animalia", "Arthropoda", "Insecta", "Diptera", "Muscidae", "Musca", "E")
  )
  standard <- c(NA, "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  ranks <- stats::setNames(rep(list(standard), length(lineages)), names(lineages))
  taxo_from_lineages(lineages, ranks = ranks)
}

# -- Ranks --------------------------------------------------------------------

test_that("the structured parser reads the rank of every node of the fixtures", {
  expected <- jsonlite::read_json(fixture("expected_ranks.json"))
  for (id in names(expected)) {
    page <- xml2::read_html(fixture(paste0("tree_", id, ".html")), encoding = "UTF-8")
    parsed <- taxodist:::.lineage_from_page(page, id)
    want <- vapply(expected[[id]], function(x) if (is.null(x)) NA_character_ else x, character(1))
    expect_identical(attr(parsed, "ranks"), want, info = paste("fixture", id))
  }

  passer <- taxodist:::.lineage_from_page(
    xml2::read_html(fixture("tree_59107.html"), encoding = "UTF-8"), "59107"
  )
  ranks <- stats::setNames(attr(passer, "ranks"), passer)
  expect_identical(
    unname(ranks[c("Animalia", "Chordata", "Reptilia", "Aves", "Passeriformes",
                   "Passeridae", "Passer", "Neornithes", "Biota")]),
    c("Kingdom", "Phylum", "Class", "Subclass", "Order", "Family", "Genus",
      "Clade", "Clade")
  )
})

test_that("ranks are inferred from the form of unlabelled names", {
  expect_identical(taxodist:::.rank_from_name("Homo sapiens"), "Species")
  expect_identical(taxodist:::.rank_from_name("Homo sapiens sapiens"), "Subspecies")
  expect_identical(taxodist:::.rank_from_name("Microtus (Agricola)"), "Subgenus")
  expect_identical(taxodist:::.rank_from_name("Natura"), NA_character_)
  expect_identical(taxodist:::.normalize_rank("[plesion] Family †", "X"), "Family")
  expect_identical(taxodist:::.normalize_rank("Taxon: Genus †", "Yi"), "Genus")
  expect_identical(taxodist:::.normalize_rank("Clade \"", "X"), "Clade")
})

test_that("the text fallback parser also records ranks", {
  fake_html <- paste0(
    '<html><body><div id="ctl00_divSubject"><b>Homo sapiens</b></div>',
    '<div id="divPageContent">Natura - nature\nClade Biota Wagner 2004\n',
    'Kingdom Animalia Linnaeus, 1758\n[crown] Clade Theria\n',
    'Genus Homo Linnaeus, 1758\nHomo sapiens Linnaeus, 1758\n</div></body></html>'
  )
  parsed <- taxodist:::.lineage_from_page(xml2::read_html(fake_html), "1")
  expect_identical(
    as.character(parsed),
    c("Natura", "Biota", "Animalia", "Theria", "Homo", "Homo sapiens")
  )
  expect_identical(
    attr(parsed, "ranks"),
    c(NA, "Clade", "Kingdom", "Clade", "Genus", "Species")
  )
})

test_that("lineage_ranks returns nodes aligned with get_lineage", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  cache <- taxodist:::.taxodist_cache
  assign("id_Passer", "59107", envir = cache)
  assign("lin_59107", c("Biota", "Animalia", "Passeridae", "Passer", "Passer domesticus"),
         envir = cache)
  assign("ranks_59107", c("Clade", "Kingdom", "Family", "Genus", "Species"), envir = cache)

  out <- lineage_ranks("Passer")
  expect_identical(out$node, get_lineage("Passer"))
  expect_identical(out$rank, c("Clade", "Kingdom", "Family", "Genus"))
  expect_identical(out$depth, 1:4)

  # A queried name absent from the page is appended with its implied rank.
  assign("id_Passer montanus", "59107", envir = cache)
  out <- lineage_ranks("Passer montanus")
  expect_identical(utils::tail(out$node, 1L), "Passer montanus")
  expect_identical(utils::tail(out$rank, 1L), "Species")
  expect_null(lineage_ranks(NA))
})

test_that("finalize keeps ranks aligned with the truncated lineage", {
  lineage <- c("Biota", "[unranked] x", "Microtus", "Microtus (Microtus)", "Microtus arvalis")
  ranks <- c("Clade", NA, "Genus", "Subgenus", "Species")
  expect_identical(
    taxodist:::.finalize_ranks(ranks, lineage, "Microtus"),
    c("Clade", "Genus")
  )
  expect_identical(
    taxodist:::.finalize_ranks(ranks, lineage, "Microtus arvalis"),
    c("Clade", "Genus", "Subgenus", "Species")
  )
})

test_that("taxo_from_lineages validates and stores ranks", {
  res <- ranked_example()
  expect_true("ranks" %in% names(res))
  expect_identical(res$ranks[[1]][[2]], "Kingdom")

  plain <- taxo_from_lineages(list(A = c("Biota", "A")))
  expect_identical(plain$ranks[[1]], c(NA_character_, NA_character_))

  expect_error(
    taxo_from_lineages(list(A = c("Biota", "A")), ranks = list(A = "Kingdom")),
    "one rank"
  )
  expect_error(
    taxo_from_lineages(list(A = c("Biota", "A")), ranks = list(B = c(NA, "Genus"))),
    "every lineage name"
  )
  expect_error(taxo_from_lineages(list(A = c("Biota", "A")), ranks = "x"), "list")
})

# -- Metrics ------------------------------------------------------------------

test_that("every metric gives the documented values", {
  res <- ranked_example()
  get <- function(metric) as.matrix(distance_matrix(res, metric = metric, progress = FALSE))

  reciprocal <- get("reciprocal")
  expect_equal(reciprocal["A", "B"], 1 / 7)
  expect_equal(reciprocal["A", "E"], 1 / 2)

  path <- get("path")
  expect_equal(path["A", "B"], 2)
  expect_equal(path["A", "C"], 6)
  expect_equal(path["A", "E"], 12)

  wu <- get("wu_palmer")
  expect_equal(wu["A", "B"], 1 - 14 / 16)
  expect_equal(wu["A", "E"], 1 - 4 / 16)

  ranked <- get("ranked")
  # Clarke & Warwick steps: genus 1, family 2, order 3, class 4, phylum 5,
  # kingdom 6 (Biota has no standard rank).
  expect_equal(ranked["A", "B"], 1)
  expect_equal(ranked["A", "C"], 3)
  expect_equal(ranked["A", "D"], 4)
  expect_equal(ranked["A", "E"], 6)

  for (m in c("reciprocal", "path", "wu_palmer", "ranked")) {
    mat <- get(m)
    expect_true(all(diag(mat) == 0), info = m)
    expect_equal(mat, t(mat), info = m)
    expect_identical(attr(distance_matrix(res, metric = m, progress = FALSE), "method"), m)
  }
})

test_that("reciprocal is ultrametric and path satisfies the triangle inequality", {
  res <- ranked_example()
  rec <- as.matrix(distance_matrix(res, progress = FALSE))
  pth <- as.matrix(distance_matrix(res, metric = "path", progress = FALSE))
  triples <- utils::combn(nrow(rec), 3)
  for (k in seq_len(ncol(triples))) {
    i <- triples[1, k]; j <- triples[2, k]; l <- triples[3, k]
    expect_lte(rec[i, l], max(rec[i, j], rec[j, l]) + 1e-12)
    expect_lte(pth[i, l], pth[i, j] + pth[j, l] + 1e-12)
  }
})

test_that("metrics treat ancestors, identical lineages and disconnected roots", {
  anc <- c("Biota", "Animalia", "Chordata")
  desc <- c("Biota", "Animalia", "Chordata", "Mammalia", "Homo")
  other <- c("Root2", "X")
  cd <- taxodist:::.compute_distance
  expect_equal(cd(anc, desc, metric = "path")$distance, 2)
  expect_equal(cd(anc, desc, metric = "wu_palmer")$distance, 1 - 6 / 8)
  expect_equal(cd(desc, desc, metric = "wu_palmer")$distance, 0)
  expect_equal(cd(desc, desc, metric = "path")$distance, 0)
  expect_identical(cd(anc, other, metric = "path")$distance, Inf)
  expect_identical(cd(anc, other, metric = "reciprocal")$distance, Inf)
  expect_equal(cd(anc, other, metric = "wu_palmer")$distance, 1)
  expect_identical(cd(anc, desc, metric = "path")$metric, "path")
})

test_that("the ranked metric uses the configured rank levels", {
  res <- ranked_example()
  old <- options(taxodist.rank_levels = c("kingdom", "phylum", "class", "order"))
  on.exit(options(old), add = TRUE)
  ranked <- as.matrix(distance_matrix(res, metric = "ranked", progress = FALSE))
  # Reduced lineages: root, Kingdom, Phylum, Class, Order, terminal.
  expect_equal(ranked["A", "B"], 1)
  expect_equal(ranked["A", "C"], 1)
  expect_equal(ranked["A", "D"], 2)
})

test_that("the ranked metric keeps taxa of different kingdoms connected", {
  res <- taxo_from_lineages(
    list(
      Homo = c("Biota", "Animalia", "Chordata", "Homo"),
      Quercus = c("Biota", "Plantae", "Tracheophyta", "Quercus")
    ),
    ranks = list(
      Homo = c("Clade", "Kingdom", "Phylum", "Genus"),
      Quercus = c("Clade", "Kingdom", "Phylum", "Genus")
    )
  )
  expect_equal(c(distance_matrix(res, metric = "ranked", progress = FALSE)), 3)
})

test_that("the ranked metric requires ranks", {
  plain <- taxo_from_lineages(list(A = c("Biota", "A"), B = c("Biota", "B")))
  expect_error(distance_matrix(plain, metric = "ranked", progress = FALSE), "requires the rank")
  expect_error(distance_matrix(plain, metric = "nope", progress = FALSE), "Unknown distance metric")
})

test_that("metric is honoured by the pairwise and analysis functions", {
  res <- ranked_example()
  cross <- cross_distances(res[1:2, ], res[3:5, ], metric = "path", progress = FALSE)
  full <- as.matrix(distance_matrix(res, metric = "path", progress = FALSE))
  expect_equal(unname(cross), unname(full[1:2, 3:5]))

  cross_ranked <- cross_distances(res[1:2, ], res[3:5, ], metric = "ranked", progress = FALSE)
  full_ranked <- as.matrix(distance_matrix(res, metric = "ranked", progress = FALSE))
  expect_equal(unname(cross_ranked), unname(full_ranked[1:2, 3:5]))

  cl <- taxo_cluster(res, metric = "path", progress = FALSE)
  expect_identical(attr(cl$dist, "method"), "path")
  expect_identical(cl$hclust$dist.method, "path")
  ord <- taxo_ordinate(res, metric = "wu_palmer", progress = FALSE)
  expect_identical(attr(ord$dist, "method"), "wu_palmer")
})

test_that("taxo_distance, closest_relative and focal_distances accept a metric", {
  lineages <- list(
    Tyrannosaurus = structure(
      c("Biota", "Dinosauria", "Theropoda", "Tyrannosaurus"),
      ranks = c("Clade", "Clade", "Order", "Genus")
    ),
    Velociraptor = structure(
      c("Biota", "Dinosauria", "Theropoda", "Velociraptor"),
      ranks = c("Clade", "Clade", "Order", "Genus")
    ),
    Triceratops = structure(
      c("Biota", "Dinosauria", "Ornithischia", "Triceratops"),
      ranks = c("Clade", "Clade", "Order", "Genus")
    )
  )
  fake <- function(taxon, verbose = FALSE) as.character(lineages[[taxon]])
  fake_ranked <- function(taxon, verbose = FALSE) lineages[[taxon]]

  mockery::stub(taxo_distance, "get_lineage", fake)
  expect_equal(taxo_distance("Tyrannosaurus", "Triceratops", metric = "path")$distance, 4)

  mockery::stub(taxo_distance, ".lineage_for_metric", function(taxon, metric, fetch, verbose) {
    fake_ranked(taxon)
  })
  expect_equal(taxo_distance("Tyrannosaurus", "Velociraptor", metric = "ranked")$distance, 1)
  expect_equal(taxo_distance("Tyrannosaurus", "Triceratops", metric = "ranked")$distance, 2)

  mockery::stub(closest_relative, "get_lineage", fake)
  cr <- closest_relative("Tyrannosaurus", c("Triceratops", "Velociraptor"), metric = "path")
  expect_identical(cr$taxon, c("Velociraptor", "Triceratops"))
  expect_equal(cr$distance, c(2, 4))

  mockery::stub(focal_distances, "get_lineage", fake)
  fd <- focal_distances("Tyrannosaurus", c("Triceratops", "Velociraptor"),
                        metric = "wu_palmer", progress = FALSE)
  expect_equal(fd$distance, c(1 - 6 / 8, 1 - 4 / 8))
})

# -- Bundles (schema 1.1) -----------------------------------------------------

test_that("bundles record the metric and ranks and survive a JSON round trip", {
  res <- ranked_example()
  for (m in c("reciprocal", "path", "wu_palmer", "ranked")) {
    bundle <- taxo_bundle(res, metric = m, progress = FALSE)
    expect_identical(bundle$schema_version, "1.1")
    expect_identical(attr(bundle$matrix, "method"), m)
    file <- tempfile(fileext = ".json")
    write_taxodist_bundle(bundle, file)
    back <- read_taxodist_bundle(file)
    unlink(file)
    expect_identical(back$metric$name, bundle$metric$name, info = m)
    expect_equal(as.matrix(back$matrix), as.matrix(bundle$matrix), info = m)
    expect_identical(back$resolution$ranks[[1]], res$ranks[[1]], info = m)
  }
  ranked <- taxo_bundle(res, metric = "ranked", progress = FALSE)
  expect_identical(ranked$metric$rank_levels, taxodist:::.default_rank_levels)
})

test_that("bundle validation recomputes distances with the stored metric", {
  res <- ranked_example()
  bundle <- taxo_bundle(res, metric = "path", progress = FALSE)
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  write_taxodist_bundle(bundle, file)
  raw <- jsonlite::read_json(file)
  raw$metric$name <- "inverse_mrca_depth"
  jsonlite::write_json(raw, file, auto_unbox = TRUE, null = "null", digits = NA)
  expect_error(read_taxodist_bundle(file), "stored distances")

  raw$metric$name <- "unknown_metric"
  jsonlite::write_json(raw, file, auto_unbox = TRUE, null = "null", digits = NA)
  expect_error(read_taxodist_bundle(file), "unknown distance metric")
})

test_that("schema 1.0 bundles are still read and must use the reciprocal metric", {
  res <- taxo_from_lineages(list(A = c("Biota", "X", "A"), B = c("Biota", "Y", "B")))
  bundle <- taxo_bundle(res, progress = FALSE)
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  write_taxodist_bundle(bundle, file)
  raw <- jsonlite::read_json(file)
  raw$schema_version <- "1.0"
  raw$taxa <- lapply(raw$taxa, function(x) { x$lineage_ranks <- NULL; x })
  jsonlite::write_json(raw, file, auto_unbox = TRUE, null = "null", digits = NA)
  old <- read_taxodist_bundle(file)
  expect_identical(old$schema_version, "1.0")
  expect_false("ranks" %in% names(old$resolution))

  raw$metric$name <- "path_length"
  raw$matrix$values <- list(list(0, 4), list(4, 0))
  jsonlite::write_json(raw, file, auto_unbox = TRUE, null = "null", digits = NA)
  expect_error(read_taxodist_bundle(file), "only the reciprocal")

  raw$schema_version <- "2.0"
  jsonlite::write_json(raw, file, auto_unbox = TRUE, null = "null", digits = NA)
  expect_error(read_taxodist_bundle(file), "Unsupported")
})

test_that("bundle JSON keeps NA ranks as null", {
  res <- ranked_example()
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  write_taxodist_bundle(taxo_bundle(res, progress = FALSE), file)
  raw <- jsonlite::read_json(file)
  expect_null(raw$taxa[[1]]$lineage_ranks[[1]])
  expect_identical(raw$taxa[[1]]$lineage_ranks[[2]], "Kingdom")
})

# -- Taxonomic distinctness ----------------------------------------------------

test_that("taxo_distinctness reproduces the hand-computed indices", {
  w <- matrix(c(0, 1, 3, 1, 0, 3, 3, 3, 0), 3, dimnames = list(letters[1:3], letters[1:3]))
  d <- stats::as.dist(w)
  out <- taxo_distinctness(c(a = 2, b = 1, c = 1), d)
  # pairs: ab = 1 (x 2*1), ac = 3 (x 2*1), bc = 3 (x 1*1); n = 4
  expect_equal(out$S, 3)
  expect_equal(out$Delta, (1 * 2 + 3 * 2 + 3 * 1) / 6)
  expect_equal(out$Delta_star, (1 * 2 + 3 * 2 + 3 * 1) / 5)
  expect_equal(out$Delta_plus, 7 / 3)
  expect_equal(out$Lambda_plus, mean((c(1, 3, 3) - 7 / 3)^2))
  expect_equal(out$SDelta_plus, 7)

  presence <- taxo_distinctness(c("a", "b", "c"), d)
  expect_equal(presence$Delta, presence$Delta_plus)
  expect_equal(presence$Delta_star, presence$Delta_plus)
  expect_s3_class(presence, "taxodist_distinctness")
  expect_invisible(print(presence))
})

test_that("taxo_distinctness agrees with vegan::taxondive", {
  skip_if_not_installed("vegan")
  res <- ranked_example()
  d <- distance_matrix(res, metric = "ranked", progress = FALSE)
  comm <- rbind(
    s1 = c(A = 5, B = 3, C = 1, D = 0, E = 2),
    s2 = c(A = 1, B = 0, C = 2, D = 4, E = 0),
    s3 = c(A = 1, B = 1, C = 1, D = 1, E = 1)
  )
  ours <- taxo_distinctness(comm, d)
  theirs <- vegan::taxondive(comm, d)
  expect_equal(ours$Delta, unname(theirs$D))
  expect_equal(ours$Delta_star, unname(theirs$Dstar))
  expect_equal(ours$Delta_plus, unname(theirs$Dplus))
  expect_equal(ours$Lambda_plus, unname(theirs$Lambda))
  expect_equal(ours$S, unname(theirs$Species))
  expect_identical(attr(ours, "metric"), "ranked")
})

test_that("taxo_distinctness handles small samples, data frames and bad input", {
  d <- stats::as.dist(matrix(c(0, 2, 2, 0), 2, dimnames = list(c("a", "b"), c("a", "b"))))
  one <- taxo_distinctness(c(a = 3, b = 0), d)
  expect_equal(one$S, 1)
  expect_true(is.na(one$Delta_plus))

  df <- data.frame(a = c(1, 0), b = c(1, 1), row.names = c("x", "y"))
  out <- taxo_distinctness(df, d)
  expect_identical(out$sample, c("x", "y"))
  expect_equal(out$Delta_plus, c(2, NA))

  inf <- stats::as.dist(matrix(c(0, Inf, Inf, 0), 2, dimnames = list(c("a", "b"), c("a", "b"))))
  expect_true(is.na(taxo_distinctness(c("a", "b"), inf)$Delta_plus))

  expect_error(taxo_distinctness(c("a", "z"), d), "no distances")
  expect_error(taxo_distinctness(c(1, 2), d), "named by taxon")
  expect_error(taxo_distinctness(c(a = -1, b = 1), d), "non-negative")
  expect_error(taxo_distinctness(list(1), d), "must be")
})

test_that("taxo_distinctness computes distances when d is not supplied", {
  res <- ranked_example()
  mockery::stub(taxo_distinctness, "distance_matrix", function(taxa, verbose, progress, metric) {
    expect_identical(metric, "ranked")
    full <- as.matrix(distance_matrix(res, metric = metric, progress = FALSE))
    structure(stats::as.dist(full[taxa, taxa]), method = metric)
  })
  d <- distance_matrix(res, metric = "ranked", progress = FALSE)
  expected <- taxo_distinctness(c("A", "B", "E"), d)
  got <- taxo_distinctness(c("A", "B", "E"))
  expect_equal(got$Delta_plus, expected$Delta_plus)
})

test_that("taxo_distinctness_test compares a sample with random draws from the pool", {
  res <- ranked_example()
  pool <- distance_matrix(res, metric = "ranked", progress = FALSE)
  set.seed(42)
  test <- taxo_distinctness_test(c("A", "B", "C"), pool, n_iter = 499)
  expect_s3_class(test, "taxodist_distinctness_test")
  stats <- test$statistics
  expect_identical(stats$statistic, c("Delta_plus", "Lambda_plus"))
  expect_equal(stats$observed[[1]], mean(c(1, 3, 3)))
  # The expected Delta+ is the mean distance of the whole pool.
  full <- as.matrix(pool)
  expect_equal(stats$expected[[1]], mean(full[lower.tri(full)]))
  expect_true(all(stats$lower <= stats$upper))
  expect_true(all(stats$p_lower > 0 & stats$p_lower <= 1))
  expect_true(all(stats$p_upper > 0 & stats$p_upper <= 1))
  expect_identical(dim(test$simulated), c(499L, 2L))
  expect_identical(test$metric, "ranked")
  # The three closest taxa form the least distinct sample.
  expect_lte(stats$p_lower[[1]], 0.2)
  expect_invisible(print(test))

  set.seed(1)
  a <- taxo_distinctness_test(c("A", "E"), pool, n_iter = 50)
  set.seed(1)
  b <- taxo_distinctness_test(c("A", "E"), pool, n_iter = 50)
  expect_identical(a, b)
})

test_that("taxo_distinctness_test validates its input", {
  pool <- distance_matrix(ranked_example(), metric = "ranked", progress = FALSE)
  expect_error(taxo_distinctness_test(c("A", "Z"), pool), "master list")
  expect_error(taxo_distinctness_test("A", pool), "between 2")
  expect_error(taxo_distinctness_test(c("A", "B"), as.matrix(pool)), "dist")
  expect_error(taxo_distinctness_test(c("A", "B"), pool, n_iter = 0), "n_iter")
  expect_error(taxo_distinctness_test(1:2, pool), "character")
  inf <- stats::as.dist(matrix(c(0, Inf, Inf, 0), 2, dimnames = list(c("a", "b"), c("a", "b"))))
  expect_error(taxo_distinctness_test(c("a", "b"), inf), "finite")
})

# -- Edge cases of the 1.0 additions -------------------------------------------

test_that("the structured parser skips rows without a name", {
  hier <- function(body) xml2::read_html(paste0(
    '<html><body><div id="ctl00_cphTabContent_divHierarchy">', body,
    '</div></body></html>'
  ))
  parsed <- taxodist:::.parse_hierarchy(hier(
    '<a href="#">Clade <b>Biota</b></a><br/><a href="#"><b> </b></a><br/>\u00a0<span>Kingdom <b>Animalia</b></span>'
  ))
  expect_identical(as.character(parsed), c("Biota", "Animalia"))
  expect_identical(attr(parsed, "ranks"), c("Clade", "Kingdom"))
  expect_null(taxodist:::.parse_hierarchy(hier("<a href='#'>no bold</a>")))
})

test_that("get_lineage_by_id stores NA ranks when the parser returns none", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  f <- get_lineage_by_id
  mockery::stub(f, "httr::RETRY", function(...) structure(list(), class = "response"))
  mockery::stub(f, "httr::status_code", function(...) 200L)
  mockery::stub(f, "httr::content", function(...) "<html></html>")
  mockery::stub(f, ".lineage_from_page", function(...) c("Biota", "X"))
  expect_identical(f(42), c("Biota", "X"))
  expect_identical(get("ranks_42", envir = taxodist:::.taxodist_cache), c(NA_character_, NA_character_))
})

test_that(".ranks_by_id refetches lineages cached without ranks", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  cache <- taxodist:::.taxodist_cache
  lin <- c("Biota", "X")

  none <- taxodist:::.ranks_by_id
  mockery::stub(none, "get_lineage_by_id", function(...) NULL)
  expect_null(none("1"))

  # A lineage restored from an old cache file has no ranks: it is fetched again.
  assign("lin_7", lin, envir = cache)
  refetch <- taxodist:::.ranks_by_id
  mockery::stub(refetch, "get_lineage_by_id", function(taxon_id, clean, verbose) {
    if (!exists("lin_7", envir = cache)) {
      assign("lin_7", lin, envir = cache)
      assign("ranks_7", c("Clade", "Genus"), envir = cache)
    }
    get("lin_7", envir = cache)
  })
  expect_identical(refetch(7), c("Clade", "Genus"))

  # The new request still yields no ranks.
  rm("lin_7", "ranks_7", envir = cache)
  assign("lin_7", lin, envir = cache)
  stale <- taxodist:::.ranks_by_id
  mockery::stub(stale, "get_lineage_by_id", function(...) lin)
  expect_null(stale("7"))

  # Ranks that do not align with the lineage are reported as NA.
  assign("ranks_7", "Clade", envir = cache)
  expect_identical(stale("7"), c(NA_character_, NA_character_))

  expect_identical(taxodist:::.cached_ranks("7", "Biota"), "Clade")
  expect_identical(taxodist:::.cached_ranks(7, lin), c(NA_character_, NA_character_))
  expect_identical(taxodist:::.cached_ranks("8", lin), c(NA_character_, NA_character_))
  expect_identical(taxodist:::.format_taxon_id(59107), "59107")
})

test_that("lineage_ranks validates input and handles missing taxa", {
  expect_error(lineage_ranks(c("a", "b")), "single taxon")
  expect_null(lineage_ranks("   "))
  no_id <- lineage_ranks
  mockery::stub(no_id, "get_taxonomicon_id", function(...) NULL)
  expect_null(no_id("Nonexistent"))
  no_ranks <- lineage_ranks
  mockery::stub(no_ranks, ".ranks_by_id", function(...) NULL)
  mockery::stub(no_ranks, "get_lineage_by_id", function(...) NULL)
  expect_null(no_ranks("59107"))
})

test_that("ranked lineages are fetched with their ranks", {
  clear_cache()
  on.exit(clear_cache(), add = TRUE)
  cache <- taxodist:::.taxodist_cache
  assign("id_Alpha", "11", envir = cache)
  assign("lin_11", c("Biota", "Animalia", "Alpha"), envir = cache)
  assign("ranks_11", c("Clade", "Kingdom", "Genus"), envir = cache)
  assign("id_Beta", "12", envir = cache)
  assign("lin_12", c("Biota", "Plantae", "Beta"), envir = cache)
  assign("ranks_12", c("Clade", "Kingdom", "Genus"), envir = cache)

  d <- distance_matrix(c("Alpha", "Beta"), metric = "ranked", progress = FALSE)
  expect_equal(c(d), 2)
  expect_identical(
    taxodist:::.lineage_for_metric("Alpha", "ranked", get_lineage),
    structure(c("Biota", "Animalia", "Alpha"), ranks = c("Clade", "Kingdom", "Genus"))
  )

  missing <- taxodist:::.ranked_lineage
  mockery::stub(missing, "lineage_ranks", function(...) NULL)
  expect_null(missing("Gamma"))
})

test_that("metric and rank-level arguments are validated", {
  expect_error(taxodist:::.match_metric(1), "must be one of")
  expect_error(taxodist:::.match_metric(NA_character_), "must be one of")
  expect_identical(taxodist:::.match_metric(c("path", "ranked")), "path")

  old <- options(taxodist.rank_levels = character(0))
  on.exit(options(old), add = TRUE)
  expect_error(taxodist:::.rank_levels(), "non-empty character")
  options(taxodist.rank_levels = c(" Order ", "Genus"))
  expect_identical(taxodist:::.rank_levels(), c("order", "genus"))

  expect_identical(taxodist:::.project_lineage(character(0), character(0), "genus"), character(0))
  projected <- taxodist:::.metric_lineages(
    list(NULL, structure(c("Biota", "Felidae", "Felis"), ranks = c("Clade", "Family", "Genus"))),
    c("x", "y"), "ranked", levels = "genus"
  )
  expect_null(projected[[1]])
  expect_identical(projected[[2]], c("Biota", "Felis"))
})

test_that("bundle metric specs and ranks are validated", {
  expect_error(
    taxodist:::.metric_from_spec(list(name = "ranked_path_length")),
    "requires .*rank_levels"
  )
  bundle <- taxo_bundle(ranked_example(), metric = "ranked", progress = FALSE)
  bundle$resolution$ranks[[1]] <- "Kingdom"
  expect_error(taxodist:::.validate_taxodist_bundle(bundle), "lineage ranks do not match")
})

test_that("taxo_distinctness reports the metric and rejects malformed samples", {
  d <- stats::as.dist(matrix(c(0, 2, 2, 0), 2, dimnames = list(c("a", "b"), c("a", "b"))))
  attr(d, "method") <- "path"
  out <- taxo_distinctness(c("a", "b"), d)
  expect_identical(attr(out, "metric"), "path")
  expect_message(expect_output(print(out), "Delta_plus"), "Metric: path")

  expect_error(taxo_distinctness(c("a", NA), d), "missing or empty")
  expect_error(taxo_distinctness(c("a", " "), d), "missing or empty")
  expect_error(taxo_distinctness(matrix(1, 1, 2), d), "unique taxa")
  expect_error(taxo_distinctness(c("a", "b"), as.matrix(d)), "dist")

  unnamed <- taxo_distinctness(matrix(1, 2, 2, dimnames = list(NULL, c("a", "b"))), d)
  expect_identical(unnamed$sample, c("sample1", "sample2"))
  single <- taxo_distinctness(matrix(1, 1, 2, dimnames = list(NULL, c("a", "b"))), d)
  expect_identical(single$sample, "sample")
})
