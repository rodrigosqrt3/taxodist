# Batch resolution, offline lineages, and portable bundles

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

test_that("print.taxodist_resolution summarizes every status", {
  resolution <- taxo_from_lineages(list(
    Alpha = c("Biota", "Animalia", "Alpha")
  ))

  expect_invisible(print(resolution))
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

