#' @importFrom httr RETRY content status_code add_headers
#' @importFrom rvest read_html html_elements html_element html_text html_attr
#' @importFrom stringr str_remove str_remove_all str_extract str_trim
NULL

# Escape every regular-expression metacharacter so that `x` is matched
# literally. The bracket expression lists `]` first so it is taken literally.
regexEscape <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)

# PCRE pattern matching `x` as a whole word. Unicode-aware lookarounds are used
# instead of `\b`, which fails next to non-word characters (e.g. a name ending
# in ")") and ignores non-ASCII letters. Use with `perl = TRUE`.
.word_pattern <- function(x) {
  paste0("(*UCP)(?<!\\w)", regexEscape(x), "(?!\\w)")
}

.contains_word <- function(text, word, ignore.case = FALSE) {
  grepl(.word_pattern(word), text, perl = TRUE, ignore.case = ignore.case)
}

# Minimum interval, in seconds, between consecutive requests to The
# Taxonomicon, so that large batches do not overload the server.
.last_request <- new.env(parent = emptyenv())

.throttle <- function() {
  interval <- getOption("taxodist.min_interval", 0.25)
  if (!is.numeric(interval) || length(interval) != 1L || is.na(interval) ||
      interval <= 0) {
    return(invisible(NULL))
  }
  last <- .last_request$time
  if (!is.null(last)) {
    wait <- interval - as.numeric(difftime(Sys.time(), last, units = "secs"))
    if (wait > 0) Sys.sleep(wait)
  }
  .last_request$time <- Sys.time()
  invisible(NULL)
}

.taxodist_user_agent <- function() {
  version <- tryCatch(
    as.character(utils::packageVersion("taxodist")),
    error = function(e) "development"
  )
  paste("taxodist R package", version)
}

# Number of attempts made for each HTTP request. Transient failures (network
# errors and 5xx responses) are retried with exponential backoff; client
# errors (4xx) are not.
.max_tries <- function() {
  tries <- getOption("taxodist.max_tries", 3L)
  if (!is.numeric(tries) || length(tries) != 1L || is.na(tries) || tries < 1) {
    return(1L)
  }
  as.integer(tries)
}

# Collapse internal runs of whitespace and trim the ends, so that
# "Homo  sapiens " and "Homo sapiens" query and cache the same name.
.normalize_taxon <- function(taxon) {
  gsub("\\s+", " ", trimws(taxon))
}

# Text scraped from The Taxonomicon is interpolated into cli messages, where
# braces would otherwise be evaluated as glue expressions.
.cli_escape <- function(x) {
  gsub("}", "}}", gsub("{", "{{", x, fixed = TRUE), fixed = TRUE)
}

.search_url <- function(taxon) {
  paste0(
    "http://taxonomicon.taxonomy.nl/TaxonList.aspx",
    "?subject=Entity&by=ScientificName&search=",
    utils::URLencode(taxon, reserved = TRUE)
  )
}

.astronomical_pattern <-
  "\\bastronomical\\b|\\bplanet\\b|\\bMinor planet\\b|\\bcomet\\b|\\basteroid\\b"

# Parse the rows of a TaxonList search page into a data frame of unique
# candidate IDs and their descriptions. Astronomical homonyms, rows without a
# valid TaxonTree link, and links without a numeric ID are skipped.
.parse_search_rows <- function(page) {
  rows <- rvest::html_elements(page, "tr")
  ids <- character(0)
  names <- character(0)

  for (row in rows) {
    text <- rvest::html_text(row, trim = TRUE)
    if (grepl(.astronomical_pattern, text, ignore.case = TRUE)) next

    links <- rvest::html_elements(row, "a[href*='TaxonTree']")
    if (length(links) == 0) next

    classes <- rvest::html_attr(links, "class")
    valid_idx <- which(!is.na(classes) & classes == "Valid")
    if (length(valid_idx) == 0) next

    id <- stringr::str_remove(
      stringr::str_extract(
        rvest::html_attr(links[[valid_idx[1]]], "href"), "id=([0-9]+)"
      ), "id="
    )
    if (is.na(id)) next

    text_entry <- trimws(gsub("\\s+", " ", text))
    text_entry <- stringr::str_remove(text_entry, "^N\\|T\\|P\\|R\\|B\\|L\\s*")

    ids <- c(ids, id)
    names <- c(names, text_entry)
  }

  keep <- !duplicated(ids)
  data.frame(id = ids[keep], name = names[keep], stringsAsFactors = FALSE)
}

# Truncate a name-resolved lineage at the queried taxon. Genus-level and
# higher queries drop nodes containing spaces (species and bracketed rank
# labels); binomial queries keep only the queried binomial. The taxon is
# matched exactly first and case-insensitively as a fallback, and appended
# only when it is absent from the page (e.g. genera whose page links only to
# their species).
.finalize_lineage <- function(lineage, taxon) {
  plan <- .finalize_plan(lineage, taxon)
  out <- lineage[plan$keep]
  if (plan$append) c(out, taxon) else out
}

# Apply the same truncation to a rank vector aligned with `lineage`. An
# appended taxon receives the rank implied by the form of its name.
.finalize_ranks <- function(ranks, lineage, taxon) {
  plan <- .finalize_plan(lineage, taxon)
  out <- ranks[plan$keep]
  if (plan$append) c(out, .rank_from_name(taxon)) else out
}

# Indices of `lineage` retained by .finalize_lineage(), and whether the
# queried taxon must be appended.
.finalize_plan <- function(lineage, taxon) {
  candidates <- which(!grepl("^\\[", lineage))
  candidates <- candidates[.keep_lineage_node(lineage[candidates], taxon)]
  kept <- lineage[candidates]
  target <- which(kept == taxon)
  if (length(target) == 0L) target <- which(tolower(kept) == tolower(taxon))
  if (length(target) > 0L) {
    list(keep = candidates[seq_len(target[[length(target)]])], append = FALSE)
  } else {
    list(keep = candidates, append = TRUE)
  }
}

# Nodes without spaces are always kept. A node with spaces (a species,
# subspecies, or subgenus name) is kept only when it is the queried taxon or
# one of its name-ancestors: the species of a queried subspecies, or the
# subgenus of a queried species. Other multi-word nodes are descendants of the
# queried taxon or residual text and are dropped.
.keep_lineage_node <- function(lineage, taxon) {
  lower <- tolower(lineage)
  target <- tolower(taxon)
  plain <- !grepl(" ", lineage)
  if (!grepl(" ", taxon)) return(plain)

  genus <- sub(" .*$", "", target)
  is_target <- lower == target
  species_of_target <- startsWith(target, paste0(lower, " "))
  subgenus_of_target <- grepl("^\\S+ \\([^()]+\\)$", lineage) &
    sub(" .*$", "", lower) == genus
  plain | is_target | species_of_target | subgenus_of_target
}

# Parse the classification tree of a TaxonTree page. Each row of
# `#ctl00_cphTabContent_divHierarchy` is an <a> (an ancestor) or a <span>
# (the page's own taxon), followed by a <br>. The taxon name is the text of
# the row's <b> element; ranks, dagger and type markers, and author strings
# live outside it, so no regular-expression cleaning is needed. Rows after the
# page's own taxon are its children and are excluded. Returns NULL when the
# page has no such hierarchy block, so the caller can fall back to text
# parsing.
.parse_hierarchy <- function(page) {
  container <- rvest::html_element(page, "#ctl00_cphTabContent_divHierarchy")
  if (inherits(container, "xml_missing") || length(container) == 0L) {
    return(NULL)
  }

  nodes <- xml2::xml_find_all(
    container,
    paste0(
      ".//node()[not(ancestor::a) and not(ancestor::span)]",
      "[self::text() or self::a or self::span or self::br]"
    )
  )

  names <- character(0)
  ranks <- character(0)
  indents <- integer(0)
  is_self <- logical(0)
  indent <- 0L
  for (node in nodes) {
    tag <- xml2::xml_name(node)
    if (tag == "text") {
      txt <- xml2::xml_text(node)
      indent <- indent + lengths(regmatches(txt, gregexpr("\u00a0", txt)))
      next
    }
    if (tag == "br") {
      indent <- 0L
      next
    }
    cls <- xml2::xml_attr(node, "class")
    if (!is.na(cls) && cls %in% c("Ordinal", "Small")) next

    bold <- xml2::xml_find_first(node, ".//b")
    if (inherits(bold, "xml_missing")) next
    name <- gsub("\\s+", " ", trimws(xml2::xml_text(bold)))
    if (!nzchar(name)) next

    names <- c(names, name)
    ranks <- c(ranks, .rank_before_name(xml2::xml_text(node), name))
    indents <- c(indents, indent)
    is_self <- c(is_self, tag == "span")
  }
  if (length(names) == 0L) return(NULL)

  self_idx <- which(is_self)
  if (length(self_idx) > 0L) {
    keep <- self_idx[[1L]]
    self_missing <- FALSE
  } else {
    # Some pages (e.g. unverified taxa) omit the row of the page's own taxon.
    # Keep the leading run of strictly increasing indentation, i.e. the path
    # from the root, and let the caller append the subject name.
    keep <- 1L
    while (keep < length(indents) && indents[[keep + 1L]] > indents[[keep]]) {
      keep <- keep + 1L
    }
    self_missing <- TRUE
  }
  out <- names[seq_len(keep)]
  attr(out, "ranks") <- ranks[seq_len(keep)]
  if (self_missing) attr(out, "self_missing") <- TRUE
  out
}

# -- Ranks --------------------------------------------------------------------

# Rank of a tree row: the text that precedes the taxon name (e.g.
# "[crown] Clade", "Superfamily \u2021?", "Anamorph genus"), without crown,
# plesion, or stem qualifiers, extinction and type markers, quotes, or a
# "Taxon:" caption. Rows without such text are species, subspecies, or
# subgenera, whose rank is implied by the form of the name.
.rank_before_name <- function(row_text, name) {
  row_text <- gsub("\\s+", " ", row_text)
  pos <- regexpr(name, row_text, fixed = TRUE)
  prefix <- if (pos > 0L) substr(row_text, 1L, pos - 1L) else ""
  .normalize_rank(prefix, name)
}

.normalize_rank <- function(prefix, name) {
  prefix <- sub("^.*:", "", prefix)
  prefix <- gsub("\\[(crown|plesion|stem)\\]", "", prefix)
  prefix <- gsub("[\u2020\u2021\u1D40?\"]", "", prefix)
  prefix <- trimws(gsub("\\s+", " ", prefix))
  if (nzchar(prefix)) prefix else .rank_from_name(name)
}

# Rank implied by the form of a scientific name: "Genus (Subgenus)" is a
# subgenus, a binomial a species, and a trinomial a subspecies. Uninomials
# carry no rank information.
.rank_from_name <- function(name) {
  if (grepl("^\\S+ \\([^()]+\\)$", name)) return("Subgenus")
  words <- length(strsplit(trimws(name), "\\s+")[[1L]])
  if (words == 2L) return("Species")
  if (words == 3L) return("Subspecies")
  NA_character_
}

# -- Internal cache -----------------------------------------------------------

.taxodist_cache <- new.env(parent = emptyenv())

.default_cache_file <- function() {
  file.path(tools::R_user_dir("taxodist", which = "cache"), "lineage-cache.rds")
}

#' Clear the taxodist lineage cache
#'
#' Clears all cached lineages stored in the current R session. Useful when
#' you suspect cached data is stale or want to force fresh retrieval.
#'
#' @return Invisibly returns `NULL`.
#' @export
#' @examples
#' \donttest{
#' clear_cache()
#' }
clear_cache <- function() {
  rm(
    list = ls(envir = .taxodist_cache, all.names = TRUE),
    envir = .taxodist_cache
  )
  invisible(NULL)
}

#' Save the taxodist lineage cache to disk
#'
#' Serialises the current session cache to an `.rds` file so it can be
#' restored in a future session with [load_cache()]. Useful for
#' reproducibility and for avoiding repeated network requests.
#'
#' @param file Path to the `.rds` file to write. If `NULL` (default), the
#'   cache is written to a persistent per-user location given by
#'   [tools::R_user_dir()], which is created if necessary.
#'
#' @return Invisibly returns the path of the written file.
#' @seealso [load_cache()], [clear_cache()]
#' @export
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".rds")
#' save_cache(tmp)
#' unlink(tmp)
#' }
save_cache <- function(file = NULL) {
  if (is.null(file)) {
    file <- .default_cache_file()
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  }
  data <- as.list(.taxodist_cache)
  n_entries <- length(data)
  data[[.parser_version_key]] <- .parser_version
  saveRDS(data, file = file)
  cli::cli_alert_success("Cache saved to {.file {file}} ({n_entries} entries).")
  invisible(file)
}

#' Load a previously saved taxodist cache from disk
#'
#' Restores lineage data saved with [save_cache()] into the current session
#' cache, avoiding network requests for taxa already retrieved in a previous
#' session. The file structure is validated before the current cache is
#' modified.
#'
#' @param file Path to an `.rds` file created by [save_cache()]. If `NULL`
#'   (default), the persistent per-user cache written by `save_cache()` is
#'   used.
#'
#' @return Invisibly returns `NULL`.
#' @seealso [save_cache()], [clear_cache()]
#' @export
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".rds")
#' save_cache(tmp)
#' load_cache(tmp)
#' unlink(tmp)
#' }
load_cache <- function(file = NULL) {
  if (is.null(file)) file <- .default_cache_file()
  if (!file.exists(file)) {
    cli::cli_abort("Cache file not found: {.file {file}}")
  }
  data <- readRDS(file)

  if (!is.list(data)) {
    cli::cli_abort("Invalid cache file: expected a named list created by {.fn save_cache}.")
  }

  data_names <- names(data)
  invalid_names <- length(data) > 0L && (
    is.null(data_names) || anyNA(data_names) ||
      any(!nzchar(data_names)) || anyDuplicated(data_names) > 0L
  )
  if (invalid_names) {
    cli::cli_abort("Invalid cache file: cache entries must have unique, non-empty names.")
  }

  # An empty cache is stored as list() and legitimately has NULL names.
  # Normalize it so the prefix checks below receive a character vector.
  if (length(data) == 0L) {
    data_names <- character(0)
  }

  id_entries <- startsWith(data_names, "id_")
  invalid_ids <- any(!vapply(data[id_entries], function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }, logical(1)))
  if (invalid_ids) {
    cli::cli_abort("Invalid cache file: taxon ID entries must be non-empty character scalars.")
  }

  lineage_entries <- startsWith(data_names, "lin_") |
    startsWith(data_names, "rawlin_")
  invalid_lineages <- any(!vapply(data[lineage_entries], is.character, logical(1)))
  if (invalid_lineages) {
    cli::cli_abort("Invalid cache file: lineage entries must be character vectors.")
  }

  matrix_lineage_entries <- startsWith(data_names, "matrix_lineage_")
  invalid_matrix_lineages <- any(!vapply(
    data[matrix_lineage_entries],
    is.character,
    logical(1)
  ))
  if (invalid_matrix_lineages) {
    cli::cli_abort(
      "Invalid cache file: matrix lineage entries must be character vectors."
    )
  }

  # Lineages saved by an older parser may contain malformed nodes. Taxon IDs
  # remain valid, so only those are restored from such a cache.
  parser_version <- data[[.parser_version_key]]
  data[[.parser_version_key]] <- NULL
  current <- is.numeric(parser_version) && length(parser_version) == 1L &&
    !is.na(parser_version) && parser_version >= .parser_version
  if (!current) {
    entry_names <- names(data)
    if (is.null(entry_names)) entry_names <- character(0)
    keep <- startsWith(entry_names, "id_")
    n_dropped <- sum(!keep)
    data <- data[keep]
    if (n_dropped > 0L) {
      cli::cli_warn(c(
        "!" = "The cache was written by an older version of taxodist.",
        "i" = paste(
          "{n_dropped} lineage entr{?y was/ies were} discarded and will be",
          "retrieved again; taxon IDs were kept."
        )
      ))
    }
  }

  list2env(data, envir = .taxodist_cache)
  cli::cli_alert_success("Cache loaded from {.file {file}} ({length(data)} entries).")
  invisible(NULL)
}

# Order of homonymous candidates: entries whose search description marks the
# name as invalid or rejected (e.g. "nom. inval.", "nom. rej.") are moved
# after the others, keeping the search order within each group.
.candidate_order <- function(descriptions) {
  invalid <- grepl("\\bnom\\.\\s*(inval|rej|nud|illeg)\\.", descriptions,
                   ignore.case = TRUE)
  order(invalid)
}

# Version of the lineage parser. Saved caches record it so that lineages
# produced by an older parser are not reused.
.parser_version <- 2L
.parser_version_key <- "meta_parser_version"

#' Inspect the current taxodist lineage cache
#'
#' Reports the number of cached entries, their total memory footprint, and
#' the names of all taxa whose lineages are stored in the current session.
#' Useful for understanding what has already been retrieved before running
#' further computations.
#'
#' @return Invisibly returns a list with elements:
#' \describe{
#'   \item{`n_lineages`}{Integer. Number of cached lineages.}
#'   \item{`n_ids`}{Integer. Number of cached taxon IDs.}
#'   \item{`taxa`}{Character vector of the taxa whose lineages are cached
#'     (the terminal node of each lineage).}
#'   \item{`ids`}{Character vector of the Taxonomicon IDs of those lineages,
#'     in the same order as `taxa`.}
#'   \item{`size_bytes`}{Numeric. Total memory used by the cache.}
#' }
#'
#' @seealso [save_cache()], [load_cache()], [clear_cache()]
#' @export
#' @examples
#' \donttest{
#' get_lineage("Tyrannosaurus")
#' get_lineage("Velociraptor")
#' cache_info()
#' }
cache_info <- function() {
  keys <- ls(envir = .taxodist_cache, all.names = TRUE)

  lin_keys <- keys[startsWith(keys, "lin_")]
  id_keys  <- keys[startsWith(keys, "id_")]

  ids <- sub("^lin_", "", lin_keys)
  taxa_names <- vapply(lin_keys, function(key) {
    lineage <- get(key, envir = .taxodist_cache, inherits = FALSE)
    if (length(lineage) == 0L) NA_character_ else lineage[[length(lineage)]]
  }, character(1), USE.NAMES = FALSE)

  size_bytes <- utils::object.size(as.list(.taxodist_cache))

  cli::cli_h2("taxodist Cache")
  cli::cli_bullets(c(
    "*" = "Lineages cached : {length(lin_keys)}",
    "*" = "IDs cached      : {length(id_keys)}",
    "*" = "Memory used     : {format(size_bytes, units = 'auto')}"
  ))

  if (length(taxa_names) > 0) {
    cli::cli_text("")
    cli::cli_text("{.strong Cached taxa:}")
    cli::cli_bullets(stats::setNames(
      paste0("{.val ", .cli_escape(taxa_names), "} (ID ", ids, ")"),
      rep(" ", length(taxa_names))
    ))
  } else {
    cli::cli_text("")
    cli::cli_alert_info("No lineages cached yet.")
  }

  invisible(list(
    n_lineages = length(lin_keys),
    n_ids      = length(id_keys),
    taxa       = taxa_names,
    ids        = ids,
    size_bytes = as.numeric(size_bytes)
  ))
}

# -- ID lookup ----------------------------------------------------------------

#' Find the Taxonomicon ID for a taxon name
#'
#' Queries The Taxonomicon (taxonomy.nl) to retrieve the internal numeric
#' identifier for a given taxon name. The search filters out non-biological
#' entities such as astronomical objects that may share the same name.
#'
#' @param taxon A character string giving the taxon name to search for.
#'   Typically a genus name (e.g., `"Tyrannosaurus"`) but species and higher
#'   ranks are also supported. Leading, trailing, and repeated internal
#'   whitespace is ignored.
#' @param verbose Logical. If `TRUE`, prints status messages during retrieval.
#'   Default is `FALSE`.
#'
#' @return A character string containing the Taxonomicon numeric ID, or `NULL`
#'   if the taxon is not found or the server cannot be reached.
#'
#' @details
#' The function queries the static search endpoint at
#' `taxonomicon.taxonomy.nl/TaxonList.aspx` and parses the resulting HTML
#' to extract the taxon ID from the hierarchy link. Only candidates whose
#' lineage descends from `Biota` are retained, which excludes astronomical
#' and other non-biological homonyms. When several biological candidates
#' remain, those whose lineage contains the queried name are preferred and a
#' warning lists the alternatives.
#'
#' Transient network failures are retried (three attempts by default; see
#' the `taxodist.max_tries` option in [taxodist-package]).
#'
#' @seealso [get_lineage()], [taxo_resolve()], [check_coverage()]
#'
#' @export
#' @examples
#' \donttest{
#' get_taxonomicon_id("Tyrannosaurus")   # returns "50841"
#' get_taxonomicon_id("Homo")
#' get_taxonomicon_id("Quercus")
#' }
get_taxonomicon_id <- function(taxon, verbose = FALSE) {
  .lookup_taxonomicon_id(taxon, verbose = verbose)$id
}

# Workhorse for get_taxonomicon_id(). Returns list(id, status) where status is
# "found", "not_found", or "retrieval_error", so that callers such as
# check_coverage() can distinguish an absent name from an unreachable server.
.lookup_taxonomicon_id <- function(taxon, verbose = FALSE) {
  taxon <- .normalize_taxon(taxon)
  cache_key <- paste0("id_", taxon)
  if (exists(cache_key, envir = .taxodist_cache, inherits = FALSE)) {
    if (verbose) cli::cli_alert_info("Using cached ID for {taxon}")
    return(list(
      id = get(cache_key, envir = .taxodist_cache, inherits = FALSE),
      status = "found"
    ))
  }

  if (verbose) cli::cli_alert_info("Searching Taxonomicon for {taxon}...")

  .throttle()
  res <- tryCatch(
    httr::RETRY(
      "GET", .search_url(taxon),
      httr::add_headers("User-Agent" = .taxodist_user_agent()),
      httr::timeout(30),
      times = .max_tries(), pause_base = 1, pause_cap = 5,
      terminate_on = 400:499, quiet = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(res) || httr::status_code(res) != 200) {
    cli::cli_warn(c(
      "!" = "Cannot reach The Taxonomicon server.",
      "x" = "The website (taxonomy.nl) appears to be offline or unreachable.",
      "i" = "Please try again later."
    ))
    return(list(id = NULL, status = "retrieval_error"))
  }

  page <- tryCatch(
    rvest::read_html(httr::content(res, "text", encoding = "UTF-8")),
    error = function(e) NULL
  )
  if (is.null(page)) {
    cli::cli_warn("Could not parse the response from The Taxonomicon.")
    return(list(id = NULL, status = "retrieval_error"))
  }

  candidates <- .parse_search_rows(page)
  bio_ids <- character(0)
  bio_text <- character(0)
  bio_lineages <- list()
  lineage_failed <- FALSE

  for (i in seq_len(nrow(candidates))) {
    candidate_lin <- get_lineage_by_id(candidates$id[[i]], clean = TRUE, verbose = FALSE)
    if (is.null(candidate_lin)) {
      lineage_failed <- TRUE
      next
    }
    if (!"Biota" %in% candidate_lin) next
    bio_ids <- c(bio_ids, candidates$id[[i]])
    bio_text <- c(bio_text, candidates$name[[i]])
    bio_lineages <- c(bio_lineages, list(candidate_lin))
  }

  if (length(bio_ids) > 1L) {
    matched <- vapply(bio_lineages, function(lin) {
      any(.contains_word(lin, taxon, ignore.case = TRUE))
    }, logical(1))
    if (any(matched)) {
      bio_ids <- bio_ids[matched]
      bio_text <- bio_text[matched]
    }
    preferred <- .candidate_order(bio_text)
    bio_ids <- bio_ids[preferred]
    bio_text <- bio_text[preferred]
  }

  if (length(bio_ids) == 0L) {
    if (lineage_failed) {
      return(list(id = NULL, status = "retrieval_error"))
    }
    if (verbose) cli::cli_alert_warning("{taxon} not found in Taxonomicon")
    return(list(id = NULL, status = "not_found"))
  }

  if (length(bio_ids) > 1L) {
    escaped <- .cli_escape(bio_text)
    warn_msg <- c(
      "!" = "Multiple valid biological entries found for {.val {taxon}}.",
      "i" = paste0("Using: ", escaped[[1]], " (ID: ", bio_ids[[1]], ")"),
      "i" = paste0(
        "To use a different entry, pass its numeric ID directly, e.g. ",
        "`get_lineage(\"", bio_ids[[2]], "\")`."
      ),
      "i" = "Other available IDs:"
    )
    for (i in 2:length(bio_ids)) {
      warn_msg <- c(warn_msg, "*" = paste0("ID ", bio_ids[[i]], ": ", escaped[[i]]))
    }
    cli::cli_warn(warn_msg)
  }

  id <- bio_ids[[1]]
  assign(cache_key, id, envir = .taxodist_cache)
  if (verbose) cli::cli_alert_success("Found {taxon} with ID {id}")
  list(id = id, status = "found")
}

# -- Lineage retrieval --------------------------------------------------------

#' Retrieve the full taxonomic lineage of a taxon
#'
#' Given a Taxonomicon numeric ID, retrieves and parses the complete
#' hierarchical lineage from root (Natura) to the taxon itself. The lineage
#' is returned as a character vector ordered from root to tip.
#'
#' @param taxon_id A single numeric or character Taxonomicon ID.
#'   Obtain this with [get_taxonomicon_id()].
#' @param clean Logical. If `TRUE` (default), removes philosophical root nodes
#'   above `Biota` (i.e., Natura, actualia, Mundus, naturalia). Dagger and
#'   superscript markers are always stripped from names.
#' @param verbose Logical. If `TRUE`, prints status messages. Default `FALSE`.
#'
#' @return A character vector of clade names from root to tip, or `NULL` if
#'   `taxon_id` is not a numeric ID or retrieval fails.
#'
#' @details
#' Lineage data is sourced from The Taxonomicon, which is based on
#' Systema Naturae 2000 (Brands, S.J., 1989 onwards). Lineage resolution varies
#' among taxonomic groups and may change when the source classification is
#' updated. Because hierarchy distances depend on lineage depth, analyses
#' should record the retrieval date and package version.
#'
#' Cleaned and uncleaned lineages are cached separately, so the value of
#' `clean` used in an earlier call never affects later calls.
#'
#' @seealso [get_lineage()], [taxo_distance()]
#'
#' @export
#' @examples
#' \donttest{
#' id <- get_taxonomicon_id("Tyrannosaurus")
#' lin <- get_lineage_by_id(id)
#' print(lin)
#' }
get_lineage_by_id <- function(taxon_id, clean = TRUE, verbose = FALSE) {
  if (is.null(taxon_id) || length(taxon_id) == 0L) {
    return(NULL)
  }
  if (length(taxon_id) > 1L) {
    cli::cli_abort(
      "{.arg taxon_id} must be a single ID, not a vector of length {length(taxon_id)}."
    )
  }
  if (is.na(taxon_id)) {
    return(NULL)
  }
  taxon_id <- .format_taxon_id(taxon_id)
  if (!grepl("^[0-9]+$", taxon_id)) {
    return(NULL)
  }

  clean_key <- paste0("lin_", taxon_id)
  raw_key <- paste0("rawlin_", taxon_id)
  cache_key <- if (clean) clean_key else raw_key
  if (exists(cache_key, envir = .taxodist_cache, inherits = FALSE)) {
    if (verbose) cli::cli_alert_info("Using cached lineage for ID {taxon_id}")
    return(get(cache_key, envir = .taxodist_cache, inherits = FALSE))
  }
  url <- paste0(
    "http://taxonomicon.taxonomy.nl/TaxonTree.aspx?id=",
    taxon_id, "&src=0"
  )
  .throttle()
  res <- tryCatch(
    httr::RETRY(
      "GET", url,
      httr::add_headers("User-Agent" = .taxodist_user_agent()),
      httr::timeout(30),
      times = .max_tries(), pause_base = 1, pause_cap = 5,
      terminate_on = 400:499, quiet = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(res) || httr::status_code(res) != 200) {
    if (verbose) cli::cli_alert_warning("Could not retrieve lineage for ID {taxon_id}")
    return(NULL)
  }
  page <- tryCatch(
    rvest::read_html(httr::content(res, "text", encoding = "UTF-8")),
    error = function(e) NULL
  )
  if (is.null(page)) {
    if (verbose) {
      cli::cli_alert_warning("Could not parse lineage for ID {taxon_id}")
    }
    return(NULL)
  }

  parsed <- .lineage_from_page(page, taxon_id)
  if (length(parsed) == 0) return(NULL)
  lineage <- as.character(parsed)
  ranks <- attr(parsed, "ranks", exact = TRUE)
  if (is.null(ranks)) ranks <- rep(NA_character_, length(lineage))

  # One request serves both variants, so cache both, with their ranks.
  start <- 1L
  idx_biota <- which(lineage == "Biota")
  if (length(idx_biota) > 0) start <- idx_biota[[1L]]
  kept <- seq.int(start, length(lineage))
  cleaned <- lineage[kept]
  assign(raw_key, lineage, envir = .taxodist_cache)
  assign(clean_key, cleaned, envir = .taxodist_cache)
  assign(paste0("rawranks_", taxon_id), ranks, envir = .taxodist_cache)
  assign(paste0("ranks_", taxon_id), ranks[kept], envir = .taxodist_cache)
  if (clean) cleaned else lineage
}

# Ranks aligned with get_lineage_by_id(taxon_id, clean). Lineages cached
# without ranks (e.g. restored from a cache file of an older version) are
# retrieved again.
.ranks_by_id <- function(taxon_id, clean = TRUE, verbose = FALSE) {
  lineage <- get_lineage_by_id(taxon_id, clean = clean, verbose = verbose)
  if (is.null(lineage)) return(NULL)
  taxon_id <- .format_taxon_id(taxon_id)
  key <- paste0(if (clean) "ranks_" else "rawranks_", taxon_id)
  if (!exists(key, envir = .taxodist_cache, inherits = FALSE)) {
    rm(
      list = intersect(
        paste0(c("lin_", "rawlin_"), taxon_id),
        ls(envir = .taxodist_cache, all.names = TRUE)
      ),
      envir = .taxodist_cache
    )
    lineage <- get_lineage_by_id(taxon_id, clean = clean, verbose = verbose)
    if (is.null(lineage) || !exists(key, envir = .taxodist_cache, inherits = FALSE)) {
      return(NULL)
    }
  }
  ranks <- get(key, envir = .taxodist_cache, inherits = FALSE)
  if (length(ranks) != length(lineage)) return(rep(NA_character_, length(lineage)))
  ranks
}

# Cached ranks of a cleaned lineage, without any network request. Returns
# NA ranks when none are cached or they do not align with `lineage`.
.cached_ranks <- function(taxon_id, lineage) {
  key <- paste0("ranks_", .format_taxon_id(taxon_id))
  if (exists(key, envir = .taxodist_cache, inherits = FALSE)) {
    ranks <- get(key, envir = .taxodist_cache, inherits = FALSE)
    if (length(ranks) == length(lineage)) return(ranks)
  }
  rep(NA_character_, length(lineage))
}

.format_taxon_id <- function(taxon_id) {
  if (is.numeric(taxon_id)) {
    taxon_id <- format(taxon_id, scientific = FALSE, trim = TRUE)
  }
  trimws(as.character(taxon_id))
}

# Extract the raw lineage (from the root of the page, before any trimming at
# Biota) of a parsed TaxonTree page. Uses the structured hierarchy block when
# present and falls back to text parsing otherwise. Repeated names are
# collapsed to their first occurrence.
.lineage_from_page <- function(page, taxon_id) {
  current_name <- rvest::html_text(
    rvest::html_element(page, "#ctl00_divSubject b"), trim = TRUE
  )
  if (!is.na(current_name)) {
    current_name <- trimws(gsub("[\u2020\u2021\u1D40?]", "", current_name))
  }

  structured <- .parse_hierarchy(page)
  if (!is.null(structured)) {
    nodes <- as.character(structured)
    ranks <- attr(structured, "ranks", exact = TRUE)
    if (isTRUE(attr(structured, "self_missing")) && !is.na(current_name) &&
        nzchar(current_name) && !identical(nodes[[length(nodes)]], current_name)) {
      subject_text <- rvest::html_text(rvest::html_element(page, "#ctl00_divSubject"))
      nodes <- c(nodes, current_name)
      ranks <- c(ranks, .rank_before_name(subject_text, current_name))
    }
  } else {
    nodes <- .parse_lineage_text(page, current_name, taxon_id)
    ranks <- attr(nodes, "ranks", exact = TRUE)
    nodes <- as.character(nodes)
  }
  first <- !duplicated(nodes)
  out <- nodes[first]
  attr(out, "ranks") <- unname(ranks[first])
  out
}

# Fallback parser for pages without the structured hierarchy block: the tree
# is read as plain text, one node per line, and rank prefixes, markers, and
# author strings are removed with regular expressions.
.parse_lineage_text <- function(page, current_name, taxon_id) {
  content_node <- rvest::html_element(page, "#divPageContent")

  if (!is.na(current_name) && !inherits(content_node, "xml_missing")) {
    tree_text <- rvest::html_text(content_node, trim = TRUE)
    raw_lines <- strsplit(tree_text, "\n")[[1]]
    raw_lines <- trimws(raw_lines)
    raw_lines <- raw_lines[nzchar(raw_lines)]
    tree_start <- grep("^Natura", raw_lines)[1]
    if (!is.na(tree_start)) raw_lines <- raw_lines[tree_start:length(raw_lines)]
    raw_lines_search <- stringr::str_remove_all(raw_lines, "[\u2020\u1D40]")
    cutoff <- which(.contains_word(raw_lines_search, current_name))[1]
    texts <- if (!is.na(cutoff)) raw_lines[seq_len(cutoff)] else raw_lines
  } else {
    links <- rvest::html_elements(page, "a[href*='TaxonTree']")
    hrefs_all <- rvest::html_attr(links, "href")
    links <- links[grepl("id=[0-9]", hrefs_all)]
    hrefs <- rvest::html_attr(links, "href")
    own_idx <- which(grepl(paste0("id=", taxon_id, "($|&)"), hrefs))
    if (length(own_idx) > 0) links <- links[seq_len(max(own_idx))]
    texts <- rvest::html_text(links, trim = TRUE)
  }
  marked <- texts |>
    stringr::str_remove_all("[\u2020\u2021\u1D40]") |>
    stringr::str_remove("^\\?\\s*")
  rank_match <- stringr::str_match(marked, paste0(
    "^(?:\\[(?:crown|plesion|stem)\\]\\s+)?",
    "(Anamorph\\s+\\w+|", .rank_alternatives, ")\\s+"
  ))[, 2L]
  lineage <- marked |>
    stringr::str_remove(paste0(
      "^\\[(crown|plesion|stem)\\]\\s+",
      "(?:(?:", .rank_alternatives, ")\\s+)?"
    )) |>
    stringr::str_remove(paste0("^(?:Anamorph\\s+\\w+|", .rank_alternatives, ")\\s+")) |>
    stringr::str_remove_all("[\u2020\u2021?]") |>
    stringr::str_trim() |>
    .extract_scientific_name()
  bare_ranks <- c(
    "Go to", "Superphylum", "Subfamily", "Suborder", "Epifamily",
    "Infraorder", "Superclass", "Subclass", "Superfamily",
    "Subgenus", "Section", "Division", "Candidatus", "Parvphylum",
    "Branch", "Supercohort", "Infracohort", "Subdivision", "Subsection",
    "Grade", "[unranked]", "(Supercluster)", "(Region)",
    "[crown]", ""
  )
  keep <- !is.na(lineage) & !lineage %in% bare_ranks &
    lineage != "" & !grepl("^\\s*$", lineage) &
    !grepl("^\"", lineage) & !grepl("^Population", lineage)
  keep[is.na(keep)] <- FALSE
  out <- lineage[keep]
  ranks <- rank_match[keep]
  missing_rank <- is.na(ranks)
  ranks[missing_rank] <- vapply(out[missing_rank], .rank_from_name, character(1))
  attr(out, "ranks") <- unname(ranks)
  out
}

.rank_alternatives <- paste(c(
  "Superdomain", "Domain", "Superkingdom", "Kingdom", "Subkingdom",
  "Infrakingdom", "Superphylum", "Phylum", "Subphylum", "Infraphylum",
  "Superclass", "Class", "Subclass", "Infraclass", "Supercohort", "Cohort",
  "Subcohort", "Magnorder", "Grandorder", "Superorder", "Order", "Suborder",
  "Infraorder", "Parvorder", "Legion", "Superfamily", "Family", "Subfamily",
  "Tribe", "Subtribe", "Genus", "Species", "Subspecies", "Clade",
  "Supergrade", "Grade", "Subgrade"
), collapse = "|")

# Lower-case words that start an author string rather than continue a
# scientific name (e.g. "de Blainville", "van der Hoeven", "R.Br. ex").
.author_particles <- c(
  "de", "del", "della", "des", "di", "do", "dos", "da", "das", "du", "d",
  "van", "von", "der", "den", "ter", "zu", "zur", "la", "le", "les", "ex",
  "in", "et", "al", "auct", "sensu", "nom", "non", "emend", "corrig", "ms",
  "pro", "sp", "spp", "var", "nec", "fide"
)

# Keep only the scientific name at the start of a cleaned tree line: the first
# word (a uninomial, possibly hyphenated), followed by up to two lower-case
# epithets for species and subspecies. Author names, years, vernacular names
# and annotations that follow are discarded, whatever their capitalization or
# alphabet. Quoted and bracketed lines are returned unchanged so the caller's
# filters can recognise them.
.extract_scientific_name <- function(x) {
  particles <- paste(.author_particles, collapse = "|")
  pattern <- paste0(
    "(*UCP)^[^\\s\",\\[(]+",
    "(?:\\s(?!(?:", particles, ")(?![\\w'\u2019-]))",
    "[a-z\u00e0-\u00ff][\\w-]*(?![\\w.'\u2019]))",
    "{0,2}"
  )
  plain <- !is.na(x) & !grepl("^[\"\\[(]", x)
  xs <- x[plain]
  m <- regexpr(pattern, xs, perl = TRUE)
  xs[m > 0] <- regmatches(xs, m)
  x[plain] <- xs
  x
}

#' Retrieve the full taxonomic lineage of a taxon by name
#'
#' A convenience wrapper that combines [get_taxonomicon_id()] and
#' [get_lineage_by_id()] into a single call. Given a taxon name, returns
#' its complete lineage from root to tip.
#'
#' @param taxon A single character string giving the taxon name or numeric
#'   Taxonomicon ID. Leading, trailing, and repeated internal whitespace is
#'   ignored, and the queried name is matched case-insensitively when
#'   truncating the lineage.
#' @param clean Logical. If `TRUE` (default), removes philosophical root nodes
#'   and cleans formatting markers.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A character vector of clade names ordered from root to tip, or
#'   `NULL` if the taxon cannot be found.
#'
#' @export
#' @examples
#' \donttest{
#' get_lineage("Tyrannosaurus")
#' get_lineage("Homo sapiens")
#' get_lineage("Quercus robur")
#' }
get_lineage <- function(taxon, clean = TRUE, verbose = FALSE) {
  if (length(taxon) != 1L) {
    cli::cli_abort("{.arg taxon} must be a single taxon name or ID.")
  }
  if (is.na(taxon)) return(NULL)
  taxon <- .normalize_taxon(as.character(taxon))
  if (!nzchar(taxon)) return(NULL)

  is_id <- grepl("^[0-9]+$", taxon)
  id <- if (is_id) taxon else get_taxonomicon_id(taxon, verbose = verbose)
  if (is.null(id)) return(NULL)

  lineage <- get_lineage_by_id(id, clean = clean, verbose = verbose)
  if (is.null(lineage)) return(NULL)

  if (!is_id) lineage <- .finalize_lineage(lineage, taxon)
  lineage
}

#' Retrieve the ranks of the nodes in a lineage
#'
#' Returns each node of the lineage of a taxon together with its rank in The
#' Taxonomicon classification (for example `"Kingdom"`, `"Clade"`,
#' `"Superfamily"`, or `"Genus"`). The nodes are exactly those returned by
#' [get_lineage()] for the same arguments.
#'
#' @inheritParams get_lineage
#'
#' @return A data frame with columns `node`, `rank`, and `depth`, one row per
#'   lineage node from root to tip, or `NULL` if the taxon cannot be found.
#'   `rank` is `NA` for nodes without a stated rank (such as the
#'   philosophical root nodes above Biota).
#'
#' @details
#' Ranks are read from the label that precedes each name in the
#' classification, without the `[crown]`, `[plesion]`, and `[stem]`
#' qualifiers or extinction markers. Rows without a label are species,
#' subspecies, or subgenera, whose rank is inferred from the form of the
#' name. Many nodes are unranked clades (`"Clade"`). Ranks are used by the
#' `"ranked"` distance method; see [taxo_distance()].
#'
#' @seealso [get_lineage()], [taxo_distance()]
#' @export
#' @examples
#' \donttest{
#' lineage_ranks("Passer")
#' }
lineage_ranks <- function(taxon, clean = TRUE, verbose = FALSE) {
  if (length(taxon) != 1L) {
    cli::cli_abort("{.arg taxon} must be a single taxon name or ID.")
  }
  if (is.na(taxon)) return(NULL)
  taxon <- .normalize_taxon(as.character(taxon))
  if (!nzchar(taxon)) return(NULL)

  is_id <- grepl("^[0-9]+$", taxon)
  id <- if (is_id) taxon else get_taxonomicon_id(taxon, verbose = verbose)
  if (is.null(id)) return(NULL)

  ranks <- .ranks_by_id(id, clean = clean, verbose = verbose)
  lineage <- get_lineage_by_id(id, clean = clean, verbose = verbose)
  if (is.null(ranks) || is.null(lineage)) return(NULL)
  if (!is_id) {
    ranks <- .finalize_ranks(ranks, lineage, taxon)
    lineage <- .finalize_lineage(lineage, taxon)
  }
  data.frame(
    node = lineage,
    rank = ranks,
    depth = seq_along(lineage),
    stringsAsFactors = FALSE
  )
}

#' Search The Taxonomicon for a taxon name
#'
#' Queries The Taxonomicon database and returns a data frame of all available
#' biological entries matching the search string. This is particularly useful
#' for exploring homonyms, historical ranks, or taxonomic synonyms before
#' computing distances.
#'
#' @param taxon A character string giving the taxon name to search for.
#' @param verbose Logical. If `TRUE`, prints status messages. Default `FALSE`.
#'
#' @return A data frame of class `"data.frame"` with columns:
#' \describe{
#'   \item{`id`}{Character. The numeric Taxonomicon ID.}
#'   \item{`name`}{Character. The full taxon description, including rank and author.}
#' }
#' Returns `NULL` if no matches are found.
#'
#' @seealso [get_lineage()], [taxo_resolve()]
#' @export
#' @examples
#' \donttest{
#' taxo_search("Bacteria")
#' taxo_search("Nereis")
#' taxo_search("Tyrannosaurus")
#' }
taxo_search <- function(taxon, verbose = FALSE) {
  .taxo_search_details(taxon, verbose = verbose)$results
}

.taxo_search_details <- function(taxon, verbose = FALSE) {
  taxon <- .normalize_taxon(taxon)

  if (verbose) cli::cli_alert_info("Searching Taxonomicon for {.val {taxon}}...")

  .throttle()
  res <- tryCatch(
    httr::RETRY(
      "GET", .search_url(taxon),
      httr::add_headers("User-Agent" = .taxodist_user_agent()),
      httr::timeout(30),
      times = .max_tries(), pause_base = 1, pause_cap = 5,
      terminate_on = 400:499, quiet = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(res) || httr::status_code(res) != 200) {
    if (verbose) cli::cli_alert_warning("Could not reach Taxonomicon")
    return(list(status = "retrieval_error", results = NULL))
  }

  page <- tryCatch(
    rvest::read_html(httr::content(res, "text", encoding = "UTF-8")),
    error = function(e) NULL
  )
  if (is.null(page)) {
    if (verbose) cli::cli_alert_warning("Could not parse the Taxonomicon response")
    return(list(status = "retrieval_error", results = NULL))
  }

  df <- .parse_search_rows(page)
  if (nrow(df) == 0L) {
    if (verbose) cli::cli_alert_warning("No matches found.")
    return(list(status = "not_found", results = NULL))
  }

  if (verbose) cli::cli_alert_success("Found {nrow(df)} entries.")
  list(status = "ok", results = df)
}

#' Resolve taxon names in a batch with an auditable result
#'
#' Resolves taxon names or numeric Taxonomicon IDs while preserving the
#' candidates considered, the selected identifier, the retrieved lineage, and
#' the resolution status for every input. This is intended for reproducible
#' workflows in which warnings emitted during a large matrix calculation would
#' otherwise be difficult to audit afterwards.
#'
#' @param taxa A character vector of taxon names or numeric Taxonomicon IDs.
#' @param ambiguity How to handle names with more than one valid biological
#'   candidate. `"warn"` (default) selects the first candidate and emits one
#'   aggregate warning, `"first"` selects it silently, and `"error"` aborts
#'   after all inputs have been inspected.
#' @param verbose Logical. If `TRUE`, prints retrieval messages.
#' @param progress Logical. If `TRUE`, displays a progress bar.
#' @param within Optional character string naming a clade (for example
#'   `"Animalia"` or `"Plantae"`). When supplied, only candidates whose
#'   lineage contains this clade are considered, which resolves many
#'   homonyms across kingdoms automatically. Matching is exact and
#'   case-insensitive. Names without any candidate inside the clade are
#'   reported as `"unresolved"`. Numeric IDs are not filtered.
#'
#' @return A data frame of class `taxodist_resolution` with one row per input
#'   and columns `input`, `resolved_name`, `id`, `status`, `n_candidates`,
#'   `lineage_depth`, `lineage`, and `candidates`. The last two are list
#'   columns. Status is one of `"resolved"`, `"ambiguous"`, `"unresolved"`,
#'   or `"retrieval_error"`.
#'
#' @details
#' Ambiguous inputs retain `status = "ambiguous"` even when the first candidate
#' is selected. Consequently, downstream code can distinguish an unambiguous
#' match from a match chosen according to the requested policy. Numeric IDs
#' bypass name search and are treated as single candidates when their lineage
#' can be retrieved.
#'
#' Names are matched after trimming and collapsing whitespace, and the stored
#' lineage is truncated at the queried taxon exactly as in [get_lineage()], so
#' `distance_matrix(taxo_resolve(x))` and `distance_matrix(x)` agree.
#'
#' Objects returned by `taxo_resolve()` can be passed directly to
#' [distance_matrix()], which reuses the stored lineages and preserves the
#' original input names as matrix labels.
#'
#' @seealso [taxo_search()], [distance_matrix()], [check_coverage()]
#' @export
#' @examples
#' \donttest{
#' resolved <- taxo_resolve(c("Tyrannosaurus", "Nereis", "50841"))
#' distance_matrix(resolved)
#'
#' # Restrict homonyms to animals
#' taxo_resolve("Nereis", within = "Animalia")
#' }
taxo_resolve <- function(taxa,
                         ambiguity = c("warn", "first", "error"),
                         verbose = FALSE,
                         progress = TRUE,
                         within = NULL) {
  ambiguity <- match.arg(ambiguity)

  if (!is.character(taxa)) {
    cli::cli_abort("{.arg taxa} must be a character vector.")
  }
  if (anyNA(taxa) || any(!nzchar(trimws(taxa)))) {
    cli::cli_abort("{.arg taxa} cannot contain missing or empty values.")
  }
  if (!is.null(within) && (!is.character(within) || length(within) != 1L ||
                           is.na(within) || !nzchar(trimws(within)))) {
    cli::cli_abort("{.arg within} must be NULL or one non-empty clade name.")
  }
  within_clade <- if (is.null(within)) NULL else tolower(.normalize_taxon(within))

  empty_candidates <- function() {
    data.frame(id = character(), name = character(), stringsAsFactors = FALSE)
  }

  failed <- function(taxon, status, candidates = empty_candidates()) {
    list(
      input = taxon,
      resolved_name = NA_character_,
      id = NA_character_,
      status = status,
      n_candidates = nrow(candidates),
      lineage_depth = NA_integer_,
      lineage = NULL,
      candidates = candidates,
      ranks = NULL
    )
  }

  resolve_one <- function(taxon) {
    is_id <- grepl("^[0-9]+$", taxon)

    if (is_id) {
      lineage <- get_lineage_by_id(taxon, clean = TRUE, verbose = verbose)
      if (is.null(lineage)) {
        return(failed(taxon, "retrieval_error"))
      }

      return(list(
        input = taxon,
        resolved_name = utils::tail(lineage, 1L),
        id = taxon,
        status = "resolved",
        n_candidates = 1L,
        lineage_depth = length(lineage),
        lineage = lineage,
        candidates = data.frame(
          id = taxon,
          name = utils::tail(lineage, 1L),
          stringsAsFactors = FALSE
        ),
        ranks = .cached_ranks(taxon, lineage)
      ))
    }

    search <- .taxo_search_details(taxon, verbose = verbose)
    if (identical(search$status, "retrieval_error")) {
      return(failed(taxon, "retrieval_error"))
    }
    candidates <- search$results
    if (identical(search$status, "not_found") || is.null(candidates) ||
        nrow(candidates) == 0L) {
      return(failed(taxon, "unresolved"))
    }

    searched_candidates <- candidates
    candidate_lineages <- lapply(candidates$id, function(id) {
      get_lineage_by_id(id, clean = TRUE, verbose = verbose)
    })
    retrieved <- !vapply(candidate_lineages, is.null, logical(1))
    valid <- vapply(candidate_lineages, function(lineage) {
      !is.null(lineage) && "Biota" %in% lineage
    }, logical(1))
    candidates <- candidates[valid, , drop = FALSE]
    candidate_lineages <- candidate_lineages[valid]

    if (nrow(candidates) == 0L) {
      # Only a failed lineage retrieval is a retrieval error; candidates that
      # were retrieved but are not biological (e.g. astronomical objects)
      # mean the name is simply not a taxon.
      # As for names outside `within`, an unresolved record carries no
      # candidates (bundle validation requires this).
      if (all(retrieved)) return(failed(taxon, "unresolved"))
      return(failed(taxon, "retrieval_error", searched_candidates))
    }

    if (!is.null(within_clade)) {
      inside <- vapply(candidate_lineages, function(lineage) {
        within_clade %in% tolower(lineage)
      }, logical(1))
      if (!any(inside)) {
        return(failed(taxon, "unresolved"))
      }
      candidates <- candidates[inside, , drop = FALSE]
      candidate_lineages <- candidate_lineages[inside]
    }

    # When a broad search returns several entries, prefer candidates whose
    # lineage contains the complete queried name. Retain all candidates if no
    # exact lineage match exists, mirroring get_taxonomicon_id()'s fallback.
    if (nrow(candidates) > 1L) {
      exact <- vapply(candidate_lineages, function(lineage) {
        any(.contains_word(lineage, taxon, ignore.case = TRUE))
      }, logical(1))
      if (any(exact)) {
        candidates <- candidates[exact, , drop = FALSE]
        candidate_lineages <- candidate_lineages[exact]
      }
      preferred <- .candidate_order(candidates$name)
      candidates <- candidates[preferred, , drop = FALSE]
      candidate_lineages <- candidate_lineages[preferred]
    }

    # Apply the same truncation as get_lineage() so that name-based and
    # resolution-based matrices are identical.
    selected_lineage <- .finalize_lineage(candidate_lineages[[1L]], taxon)
    selected_ranks <- .finalize_ranks(
      .cached_ranks(candidates$id[[1L]], candidate_lineages[[1L]]),
      candidate_lineages[[1L]],
      taxon
    )
    n_candidates <- nrow(candidates)
    rownames(candidates) <- NULL
    list(
      input = taxon,
      resolved_name = utils::tail(selected_lineage, 1L),
      id = candidates$id[[1L]],
      status = if (n_candidates > 1L) "ambiguous" else "resolved",
      n_candidates = n_candidates,
      lineage_depth = length(selected_lineage),
      lineage = selected_lineage,
      candidates = candidates,
      ranks = selected_ranks
    )
  }

  keys <- .normalize_taxon(taxa)
  unique_keys <- unique(keys)
  progress_id <- NULL
  if (progress && length(unique_keys) > 0L) {
    progress_id <- cli::cli_progress_bar(
      "Resolving taxa",
      total = length(unique_keys)
    )
  }
  resolved_unique <- lapply(unique_keys, function(taxon) {
    result <- resolve_one(taxon)
    if (progress) cli::cli_progress_update(id = progress_id)
    result
  })
  if (progress && length(unique_keys) > 0L) {
    cli::cli_progress_done(id = progress_id)
  }
  resolved <- resolved_unique[match(keys, unique_keys)]

  result <- data.frame(
    input = as.character(taxa),
    resolved_name = vapply(resolved, `[[`, character(1), "resolved_name"),
    id = vapply(resolved, `[[`, character(1), "id"),
    status = vapply(resolved, `[[`, character(1), "status"),
    n_candidates = vapply(resolved, `[[`, integer(1), "n_candidates"),
    lineage_depth = vapply(resolved, `[[`, integer(1), "lineage_depth"),
    stringsAsFactors = FALSE
  )
  result$lineage <- I(lapply(resolved, `[[`, "lineage"))
  result$candidates <- I(lapply(resolved, `[[`, "candidates"))
  result$ranks <- I(lapply(resolved, `[[`, "ranks"))
  class(result) <- c("taxodist_resolution", "data.frame")
  attr(result, "source") <- "The Taxonomicon"
  attr(result, "source_url") <- "http://taxonomicon.taxonomy.nl"
  attr(result, "retrieved_at") <- format(Sys.time(), tz = "UTC", usetz = TRUE)

  ambiguous_inputs <- result$input[result$status == "ambiguous"]
  if (length(ambiguous_inputs) > 0L) {
    details <- paste0(
      ambiguous_inputs,
      " (",
      result$n_candidates[result$status == "ambiguous"],
      " candidates)"
    )
    if (ambiguity == "error") {
      cli::cli_abort(c(
        "x" = "Ambiguous taxon names were found.",
        "i" = "{paste(details, collapse = ', ')}",
        "i" = "Inspect {.fn taxo_search} and pass numeric IDs to resolve them explicitly."
      ))
    }
    if (ambiguity == "warn") {
      cli::cli_warn(c(
        "!" = "Ambiguous taxon names were resolved using the first candidate.",
        "i" = "{paste(details, collapse = ', ')}",
        "i" = "Inspect the {.field candidates} column, use {.arg within}, or pass numeric IDs explicitly."
      ))
    }
  }

  result
}

#' Create an auditable resolution from user-supplied lineages
#'
#' Builds a `taxodist_resolution` object without consulting an online
#' taxonomy service. This is useful for unpublished classifications, curated
#' local taxonomies, frozen analyses, and fully offline workflows.
#'
#' @param lineages A named list of character vectors ordered from root to the
#'   focal taxon. List names become the input labels.
#' @param ids Optional character vector of unique identifiers, in the same
#'   order as `lineages`. A named vector is matched by name. If omitted,
#'   deterministic identifiers of the form `custom:<input label>` are created.
#' @param source A non-empty character string describing the lineage source.
#' @param ranks Optional named list of character vectors giving the rank of
#'   every node, aligned with `lineages` (use `NA` for unranked nodes). Ranks
#'   are required only by the `"ranked"` distance method. A list without
#'   names is matched by position.
#'
#' @return A `taxodist_resolution` object that can be passed to
#'   [distance_matrix()] or [taxo_bundle()].
#'
#' @seealso [taxo_resolve()], [distance_matrix()], [taxo_bundle()]
#' @export
#' @examples
#' lineages <- list(
#'   Alpha = c("Biota", "Animalia", "Alpha"),
#'   Beta = c("Biota", "Animalia", "Beta")
#' )
#' resolved <- taxo_from_lineages(lineages)
#' distance_matrix(resolved)
#'
#' # With ranks, the rank-based distance can be used
#' ranked <- taxo_from_lineages(
#'   lineages,
#'   ranks = list(
#'     Alpha = c(NA, "Kingdom", "Genus"),
#'     Beta = c(NA, "Kingdom", "Genus")
#'   )
#' )
#' distance_matrix(ranked, metric = "ranked")
taxo_from_lineages <- function(lineages, ids = NULL, source = "user-supplied",
                               ranks = NULL) {
  if (!is.list(lineages)) {
    cli::cli_abort("{.arg lineages} must be a named list of character vectors.")
  }
  labels <- names(lineages)
  if (is.null(labels) || length(labels) != length(lineages) ||
      anyNA(labels) || any(!nzchar(trimws(labels))) || anyDuplicated(labels)) {
    cli::cli_abort(
      "{.arg lineages} must have unique, non-empty names for every entry."
    )
  }
  valid_lineage <- vapply(lineages, function(lineage) {
    is.character(lineage) && length(lineage) > 0L && !anyNA(lineage) &&
      all(nzchar(trimws(lineage)))
  }, logical(1))
  if (any(!valid_lineage)) {
    cli::cli_abort(
      "Every lineage must be a non-empty character vector without missing or empty nodes."
    )
  }
  if (!is.character(source) || length(source) != 1L || is.na(source) ||
      !nzchar(trimws(source))) {
    cli::cli_abort("{.arg source} must be one non-empty character string.")
  }

  if (is.null(ids)) {
    ids <- paste0("custom:", labels)
  } else {
    if (!is.character(ids)) {
      cli::cli_abort("{.arg ids} must be a character vector.")
    }
    if (!is.null(names(ids))) {
      if (!all(labels %in% names(ids))) {
        cli::cli_abort("Named {.arg ids} must contain every lineage name.")
      }
      ids <- unname(ids[labels])
    }
    if (length(ids) != length(lineages) || anyNA(ids) ||
        any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
      cli::cli_abort(
        "{.arg ids} must contain one unique, non-empty identifier per lineage."
      )
    }
  }

  if (is.null(ranks)) {
    ranks <- lapply(lineages, function(lineage) rep(NA_character_, length(lineage)))
  } else {
    if (!is.list(ranks)) {
      cli::cli_abort("{.arg ranks} must be a list of character vectors.")
    }
    if (!is.null(names(ranks))) {
      if (!all(labels %in% names(ranks))) {
        cli::cli_abort("Named {.arg ranks} must contain every lineage name.")
      }
      ranks <- ranks[labels]
    }
    valid_ranks <- length(ranks) == length(lineages) &&
      all(vapply(seq_along(ranks), function(i) {
        r <- ranks[[i]]
        (is.character(r) || all(is.na(r))) && length(r) == length(lineages[[i]])
      }, logical(1)))
    if (!valid_ranks) {
      cli::cli_abort(
        "{.arg ranks} must give one rank (or {.val NA}) per node of every lineage."
      )
    }
    ranks <- lapply(ranks, as.character)
  }

  resolved_names <- vapply(lineages, utils::tail, character(1), n = 1L)
  result <- data.frame(
    input = labels,
    resolved_name = resolved_names,
    id = unname(ids),
    status = rep("resolved", length(lineages)),
    n_candidates = rep(1L, length(lineages)),
    lineage_depth = lengths(lineages),
    stringsAsFactors = FALSE
  )
  result$lineage <- I(unname(lineages))
  result$candidates <- I(lapply(seq_along(lineages), function(i) {
    data.frame(
      id = ids[[i]],
      name = resolved_names[[i]],
      stringsAsFactors = FALSE
    )
  }))
  result$ranks <- I(unname(ranks))
  class(result) <- c("taxodist_resolution", "data.frame")
  attr(result, "source") <- source
  attr(result, "source_url") <- NA_character_
  attr(result, "retrieved_at") <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  result
}
