#' Create a portable taxodist analysis bundle
#'
#' Combines an auditable taxon-resolution table, its distance matrix, the
#' distance definition, and provenance metadata in one self-contained object.
#' The bundle can be written to a language-independent JSON file with
#' [write_taxodist_bundle()].
#'
#' @param taxa A character vector accepted by [taxo_resolve()] or an existing
#'   `taxodist_resolution` object.
#' @param ambiguity,verbose,progress Passed to [taxo_resolve()] when `taxa` is
#'   a character vector.
#'
#' @return An object of class `taxodist_bundle` containing `schema_version`,
#'   `created_at`, `source`, `software`, `metric`, `resolution`, and `matrix`.
#'
#' @export
#' @examples
#' \donttest{
#' bundle <- taxo_bundle(c("Tyrannosaurus", "Triceratops", "Homo"))
#' distance_matrix(bundle)
#' }
taxo_bundle <- function(taxa,
                        ambiguity = c("warn", "first", "error"),
                        verbose = FALSE,
                        progress = TRUE) {
  ambiguity <- match.arg(ambiguity)
  resolution <- if (inherits(taxa, "taxodist_resolution")) {
    taxa
  } else {
    taxo_resolve(
      taxa,
      ambiguity = ambiguity,
      verbose = verbose,
      progress = progress
    )
  }

  matrix <- distance_matrix(resolution, progress = progress)
  created_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  package_version <- tryCatch(
    as.character(utils::packageVersion("taxodist")),
    error = function(e) NA_character_
  )
  source_name <- attr(resolution, "source", exact = TRUE)
  if (is.null(source_name) || is.na(source_name) || !nzchar(source_name)) {
    source_name <- "The Taxonomicon"
  }
  source_url <- attr(resolution, "source_url", exact = TRUE)
  if (is.null(source_url)) {
    source_url <- if (identical(source_name, "The Taxonomicon")) {
      "http://taxonomicon.taxonomy.nl"
    } else {
      NA_character_
    }
  }

  structure(
    list(
      schema_version = "1.0",
      created_at = created_at,
      source = list(
        name = source_name,
        url = source_url,
        retrieved_at = attr(resolution, "retrieved_at", exact = TRUE)
      ),
      software = list(
        name = "taxodist",
        version = package_version,
        language = "R"
      ),
      metric = list(
        name = "inverse_mrca_depth",
        definition = "0 for identical lineages; otherwise 1 / depth(MRCA)",
        root_depth = 1L,
        common_ancestor = "continuous common lineage prefix"
      ),
      resolution = resolution,
      matrix = matrix
    ),
    class = "taxodist_bundle"
  )
}

.validate_taxodist_bundle <- function(bundle) {
  if (!inherits(bundle, "taxodist_bundle") || !is.list(bundle)) {
    cli::cli_abort("{.arg bundle} must be a {.cls taxodist_bundle} object.")
  }
  required <- c(
    "schema_version", "created_at", "source", "software", "metric",
    "resolution", "matrix"
  )
  if (!all(required %in% names(bundle))) {
    cli::cli_abort("Invalid {.cls taxodist_bundle}: required fields are missing.")
  }
  if (!identical(bundle$schema_version, "1.0")) {
    cli::cli_abort(
      "Unsupported taxodist bundle schema: {.val {bundle$schema_version}}."
    )
  }
  if (!inherits(bundle$resolution, "taxodist_resolution")) {
    cli::cli_abort("Invalid bundle: {.field resolution} has the wrong class.")
  }
  resolution_fields <- c(
    "input", "resolved_name", "id", "status", "n_candidates",
    "lineage_depth", "lineage", "candidates"
  )
  if (!all(resolution_fields %in% names(bundle$resolution))) {
    cli::cli_abort("Invalid bundle: resolution fields are missing.")
  }
  allowed_status <- c(
    "resolved", "ambiguous", "unresolved", "retrieval_error"
  )
  if (any(!bundle$resolution$status %in% allowed_status)) {
    cli::cli_abort("Invalid bundle: unknown resolution status.")
  }
  if (anyNA(bundle$resolution$n_candidates) ||
      any(bundle$resolution$n_candidates < 0L)) {
    cli::cli_abort("Invalid bundle: candidate counts must be non-negative integers.")
  }
  for (i in seq_len(nrow(bundle$resolution))) {
    lineage <- bundle$resolution$lineage[[i]]
    candidates <- bundle$resolution$candidates[[i]]
    status <- bundle$resolution$status[[i]]
    if (!is.data.frame(candidates) ||
        !all(c("id", "name") %in% names(candidates))) {
      cli::cli_abort("Invalid bundle: malformed candidate table at row {i}.")
    }
    if (nrow(candidates) != bundle$resolution$n_candidates[[i]]) {
      cli::cli_abort("Invalid bundle: candidate count mismatch at row {i}.")
    }
    if (is.null(lineage)) {
      if (!is.na(bundle$resolution$lineage_depth[[i]])) {
        cli::cli_abort("Invalid bundle: lineage depth mismatch at row {i}.")
      }
    } else {
      if (!is.character(lineage) || anyNA(lineage) ||
          any(!nzchar(trimws(lineage)))) {
        cli::cli_abort("Invalid bundle: malformed lineage at row {i}.")
      }
      if (length(lineage) != bundle$resolution$lineage_depth[[i]]) {
        cli::cli_abort("Invalid bundle: lineage depth mismatch at row {i}.")
      }
    }
    if (status %in% c("resolved", "ambiguous")) {
      if (is.na(bundle$resolution$id[[i]]) || is.null(lineage) ||
          nrow(candidates) == 0L) {
        cli::cli_abort("Invalid bundle: incomplete resolved record at row {i}.")
      }
      if (!identical(as.character(candidates$id[[1L]]), bundle$resolution$id[[i]])) {
        cli::cli_abort("Invalid bundle: selected candidate mismatch at row {i}.")
      }
    } else {
      if (!is.na(bundle$resolution$id[[i]]) || !is.null(lineage) ||
          !is.na(bundle$resolution$resolved_name[[i]])) {
        cli::cli_abort("Invalid bundle: incomplete unresolved record at row {i}.")
      }
      if (identical(status, "unresolved") && nrow(candidates) != 0L) {
        cli::cli_abort("Invalid bundle: unresolved record has candidates at row {i}.")
      }
    }
  }
  if (!inherits(bundle$matrix, "dist")) {
    cli::cli_abort("Invalid bundle: {.field matrix} must be a {.cls dist} object.")
  }
  labels <- attr(bundle$matrix, "Labels", exact = TRUE)
  if (!identical(labels, bundle$resolution$input)) {
    cli::cli_abort(
      "Invalid bundle: matrix labels do not match the resolution inputs."
    )
  }
  expected_matrix <- distance_matrix(bundle$resolution, progress = FALSE)
  if (!isTRUE(all.equal(
    as.matrix(bundle$matrix),
    as.matrix(expected_matrix),
    tolerance = sqrt(.Machine$double.eps),
    check.attributes = FALSE
  ))) {
    cli::cli_abort(
      "Invalid bundle: stored distances do not match the stored lineages."
    )
  }
  invisible(bundle)
}

.bundle_taxon_records <- function(resolution) {
  lapply(seq_len(nrow(resolution)), function(i) {
    candidate_table <- resolution$candidates[[i]]
    candidates <- if (is.null(candidate_table) || nrow(candidate_table) == 0L) {
      list()
    } else {
      lapply(seq_len(nrow(candidate_table)), function(j) {
        list(
          id = candidate_table$id[[j]],
          name = candidate_table$name[[j]]
        )
      })
    }
    list(
      input = resolution$input[[i]],
      resolved_name = resolution$resolved_name[[i]],
      id = resolution$id[[i]],
      status = resolution$status[[i]],
      n_candidates = resolution$n_candidates[[i]],
      lineage_depth = resolution$lineage_depth[[i]],
      lineage = if (is.null(resolution$lineage[[i]])) {
        NULL
      } else {
        I(resolution$lineage[[i]])
      },
      candidates = candidates
    )
  })
}

.bundle_matrix_rows <- function(matrix) {
  full <- as.matrix(matrix)
  lapply(seq_len(nrow(full)), function(i) {
    lapply(as.numeric(full[i, , drop = TRUE]), function(value) {
      if (is.na(value)) return(NULL)
      if (is.infinite(value)) {
        return(if (value > 0) "Infinity" else "-Infinity")
      }
      value
    })
  })
}

#' Write a taxodist bundle as portable JSON
#'
#' @param bundle A `taxodist_bundle` created by [taxo_bundle()].
#' @param file Path to the JSON file to create.
#' @param pretty Logical. If `TRUE` (default), indent the JSON for readability.
#'
#' @return Invisibly returns the normalized output path.
#' @export
write_taxodist_bundle <- function(bundle, file, pretty = TRUE) {
  .validate_taxodist_bundle(bundle)
  full <- as.matrix(bundle$matrix)
  portable <- list(
    format = "taxodist_bundle",
    schema_version = bundle$schema_version,
    created_at = bundle$created_at,
    source = bundle$source,
    software = bundle$software,
    metric = bundle$metric,
    taxa = .bundle_taxon_records(bundle$resolution),
    matrix = list(
      labels = I(rownames(full)),
      values = .bundle_matrix_rows(bundle$matrix)
    )
  )
  jsonlite::write_json(
    portable,
    path = file,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = NA,
    pretty = pretty
  )
  invisible(normalizePath(file, winslash = "/", mustWork = TRUE))
}

.json_character <- function(value) {
  if (is.null(value) || length(value) == 0L) NA_character_ else as.character(value)
}

.json_integer <- function(value) {
  if (is.null(value) || length(value) == 0L) NA_integer_ else as.integer(value)
}

.resolution_from_json <- function(records, source, source_url, retrieved_at) {
  if (length(records) == 0L) {
    result <- data.frame(
      input = character(),
      resolved_name = character(),
      id = character(),
      status = character(),
      n_candidates = integer(),
      lineage_depth = integer(),
      stringsAsFactors = FALSE
    )
    result$lineage <- I(list())
    result$candidates <- I(list())
  } else {
    result <- data.frame(
      input = vapply(records, function(x) .json_character(x$input), character(1)),
      resolved_name = vapply(
        records,
        function(x) .json_character(x$resolved_name),
        character(1)
      ),
      id = vapply(records, function(x) .json_character(x$id), character(1)),
      status = vapply(records, function(x) .json_character(x$status), character(1)),
      n_candidates = vapply(
        records,
        function(x) .json_integer(x$n_candidates),
        integer(1)
      ),
      lineage_depth = vapply(
        records,
        function(x) .json_integer(x$lineage_depth),
        integer(1)
      ),
      stringsAsFactors = FALSE
    )
    result$lineage <- I(lapply(records, function(x) {
      if (is.null(x$lineage)) NULL else as.character(unlist(x$lineage))
    }))
    result$candidates <- I(lapply(records, function(x) {
      if (is.null(x$candidates) || length(x$candidates) == 0L) {
        return(data.frame(
          id = character(),
          name = character(),
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        id = vapply(x$candidates, function(y) .json_character(y$id), character(1)),
        name = vapply(x$candidates, function(y) .json_character(y$name), character(1)),
        stringsAsFactors = FALSE
      )
    }))
  }
  class(result) <- c("taxodist_resolution", "data.frame")
  attr(result, "source") <- source
  attr(result, "source_url") <- source_url
  attr(result, "retrieved_at") <- retrieved_at
  result
}

.matrix_from_json <- function(matrix_data) {
  labels <- as.character(unlist(matrix_data$labels))
  n <- length(labels)
  if (n == 0L) {
    return(structure(
      numeric(),
      Size = 0L,
      Labels = character(),
      Diag = FALSE,
      Upper = FALSE,
      class = "dist"
    ))
  }
  rows <- matrix_data$values
  if (length(rows) != n) {
    cli::cli_abort("Invalid bundle JSON: matrix row count does not match labels.")
  }
  values <- vapply(rows, function(row) {
    if (length(row) != n) {
      cli::cli_abort("Invalid bundle JSON: matrix must be square.")
    }
    vapply(row, function(value) {
      if (is.null(value)) return(NA_real_)
      if (identical(value, "Infinity")) return(Inf)
      if (identical(value, "-Infinity")) return(-Inf)
      as.numeric(value)
    }, numeric(1))
  }, numeric(n))
  full <- t(values)
  dimnames(full) <- list(labels, labels)
  if (anyNA(diag(full)) || any(diag(full) != 0)) {
    cli::cli_abort("Invalid bundle JSON: matrix diagonal must contain zeros.")
  }
  transpose <- t(full)
  symmetric <- (is.na(full) & is.na(transpose)) |
    (!is.na(full) & !is.na(transpose) & full == transpose)
  if (!all(symmetric)) {
    cli::cli_abort("Invalid bundle JSON: distance matrix must be symmetric.")
  }
  stats::as.dist(full)
}

#' Read a portable taxodist JSON bundle
#'
#' @param file Path to a JSON file created by [write_taxodist_bundle()].
#'
#' @return A reconstructed `taxodist_bundle` object.
#' @export
read_taxodist_bundle <- function(file) {
  if (!file.exists(file)) {
    cli::cli_abort("Bundle file not found: {.file {file}}")
  }
  raw <- tryCatch(
    jsonlite::read_json(file, simplifyVector = FALSE),
    error = function(e) {
      cli::cli_abort("Could not parse bundle JSON: {conditionMessage(e)}")
    }
  )
  if (!identical(raw$format, "taxodist_bundle")) {
    cli::cli_abort("Invalid bundle JSON: unrecognized format.")
  }
  if (!identical(raw$schema_version, "1.0")) {
    cli::cli_abort(
      "Unsupported taxodist bundle schema: {.val {raw$schema_version}}."
    )
  }
  required <- c("created_at", "source", "software", "metric", "taxa", "matrix")
  if (!all(required %in% names(raw))) {
    cli::cli_abort("Invalid bundle JSON: required fields are missing.")
  }

  source_name <- .json_character(raw$source$name)
  source_url <- .json_character(raw$source$url)
  retrieved_at <- .json_character(raw$source$retrieved_at)
  resolution <- .resolution_from_json(
    raw$taxa,
    source_name,
    source_url,
    retrieved_at
  )
  matrix <- .matrix_from_json(raw$matrix)
  bundle <- structure(
    list(
      schema_version = raw$schema_version,
      created_at = .json_character(raw$created_at),
      source = list(
        name = source_name,
        url = source_url,
        retrieved_at = retrieved_at
      ),
      software = raw$software,
      metric = raw$metric,
      resolution = resolution,
      matrix = matrix
    ),
    class = "taxodist_bundle"
  )
  .validate_taxodist_bundle(bundle)
  bundle
}

#' Print a taxodist bundle
#'
#' @param x A `taxodist_bundle` object.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns `x`.
#' @export
print.taxodist_bundle <- function(x, ...) {
  .validate_taxodist_bundle(x)
  counts <- table(factor(
    x$resolution$status,
    levels = c("resolved", "ambiguous", "unresolved", "retrieval_error")
  ))
  cli::cli_h2("taxodist Bundle")
  cli::cli_bullets(c(
    "*" = "Schema: {x$schema_version}",
    "*" = "Created: {x$created_at}",
    "*" = "Taxa: {nrow(x$resolution)}",
    "*" = "Resolved: {counts[['resolved']]}",
    "*" = "Ambiguous: {counts[['ambiguous']]}",
    "*" = "Unresolved: {counts[['unresolved']]}",
    "*" = "Retrieval errors: {counts[['retrieval_error']]}"
  ))
  invisible(x)
}
