#' @importFrom httr RETRY content status_code add_headers
#' @importFrom rvest read_html html_elements html_element html_text html_attr
#' @importFrom stringr str_remove str_remove_all str_extract str_trim
NULL

regexEscape <- function(x) gsub("([.\\^$*+?{}\\[\\]|()])", "\\\\\\1", x)

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
  match_target <- function(lineage) {
    idx <- which(lineage == taxon)
    if (length(idx) == 0L) idx <- which(tolower(lineage) == tolower(taxon))
    idx
  }

  if (!grepl(" ", taxon)) {
    lineage <- lineage[!grepl(" ", lineage)]
    lineage <- lineage[!grepl("^\\[", lineage)]
    target_idx <- match_target(lineage)
    if (length(target_idx) > 0) {
      lineage <- lineage[seq_len(target_idx[length(target_idx)])]
    } else {
      lineage <- c(lineage, taxon)
    }
  } else {
    lineage <- lineage[!grepl(" ", lineage) | tolower(lineage) == tolower(taxon)]
    target_idx <- match_target(lineage)
    if (length(target_idx) > 0) {
      lineage <- lineage[seq_len(target_idx[1])]
    } else {
      lineage <- c(lineage, taxon)
    }
  }
  lineage
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
  saveRDS(data, file = file)
  cli::cli_alert_success("Cache saved to {.file {file}} ({length(data)} entries).")
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

  list2env(data, envir = .taxodist_cache)
  cli::cli_alert_success("Cache loaded from {.file {file}} ({length(data)} entries).")
  invisible(NULL)
}

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
    pattern <- paste0("\\b", regexEscape(taxon), "\\b")
    matched <- vapply(bio_lineages, function(lin) {
      any(grepl(pattern, lin, ignore.case = TRUE))
    }, logical(1))
    if (any(matched)) {
      bio_ids <- bio_ids[matched]
      bio_text <- bio_text[matched]
    }
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
  if (is.numeric(taxon_id)) {
    taxon_id <- format(taxon_id, scientific = FALSE, trim = TRUE)
  }
  taxon_id <- trimws(as.character(taxon_id))
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

  current_name <- rvest::html_text(
    rvest::html_element(page, "#ctl00_divSubject b"), trim = TRUE
  )
  content_node <- rvest::html_element(page, "#divPageContent")

  if (!is.na(current_name) && !is.null(content_node)) {
    tree_text <- rvest::html_text(content_node, trim = TRUE)
    raw_lines <- strsplit(tree_text, "\n")[[1]]
    raw_lines <- trimws(raw_lines)
    raw_lines <- raw_lines[nzchar(raw_lines)]
    tree_start <- grep("^Natura", raw_lines)[1]
    if (!is.na(tree_start)) raw_lines <- raw_lines[tree_start:length(raw_lines)]
    raw_lines_search <- stringr::str_remove_all(raw_lines, "[\u2020\u1D40]")
    cutoff <- grep(paste0("\\b", regexEscape(current_name), "\\b"), raw_lines_search)[1]
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
  lineage <- texts |>
    stringr::str_remove_all("[\u2020\u1D40]") |>
    stringr::str_remove("^\\[crown\\]\\s+(Clade|Grandorder|Order|Superorder|Infraorder|Suborder|Class|Superclass|Subclass|Infraclass|Family|Superfamily|Subfamily|Tribe|Subtribe|Kingdom|Subkingdom|Infrakingdom|Domain|Superkingdom|Phylum|Subphylum|Genus|Species)?\\s*") |>
    stringr::str_remove("^(Clade |Kingdom |Phylum |Superphylum |Subphylum |Infraphylum |Class |Order |Suborder |Infraorder |Parvorder |Grandorder |Magnorder |Cohort |Subcohort |Legion |Family |Subfamily |Tribe |Subtribe |Genus |Species |Subkingdom |Infrakingdom |Superclass |Subclass |Infraclass |Superorder |Superfamily |Domain |Superkingdom |Grade |Subgrade |Supergrade )") |>
    stringr::str_remove("\\s+[A-Z][a-z\u00e1\u00e0\u00e2\u00e3\u00e9\u00e8\u00ea\u00ed\u00ef\u00f3\u00f4\u00f5\u00f6\u00fa\u00fc\u00e7].*$") |>
    stringr::str_remove("\\s+[A-Z]\\.[A-Z]\\..*$") |>
    stringr::str_remove("\\s+auct\\..*$") |>
    stringr::str_remove("\\s+von.*$") |>
    stringr::str_remove("\\s+\\([A-Z][a-z].*$") |>
    stringr::str_remove("\\s+\\(\\d{4}\\).*$") |>
    stringr::str_remove("\\s+\\[.*$") |>
    stringr::str_remove("\\s+[A-Z]\\.$") |>
    stringr::str_remove("\\s+\\([a-z].*$") |>
    stringr::str_remove('\\s+".*$') |>
    stringr::str_remove("^\".*") |>
    stringr::str_trim()
  bare_ranks <- c(
    "Go to", "Superphylum", "Subfamily", "Suborder", "Epifamily",
    "Infraorder", "Superclass", "Subclass", "Superfamily",
    "Subgenus", "Section", "Division", "Candidatus", "Parvphylum",
    "Branch", "Supercohort", "Infracohort", "Subdivision", "Subsection",
    "Grade", "[unranked]", "(Supercluster)", "(Region)",
    "[crown]", ""
  )
  lineage <- lineage[!lineage %in% bare_ranks]
  lineage <- lineage[lineage != "" & !grepl("^\\s*$", lineage)]
  lineage <- lineage[!grepl("^\"", lineage)]
  lineage <- lineage[!grepl("^Population", lineage)]
  lineage <- unique(lineage)
  if (length(lineage) == 0) return(NULL)

  # One request serves both variants, so cache both.
  cleaned <- lineage
  idx_biota <- which(lineage == "Biota")
  if (length(idx_biota) > 0) {
    cleaned <- lineage[idx_biota[1]:length(lineage)]
  }
  assign(raw_key, lineage, envir = .taxodist_cache)
  assign(clean_key, cleaned, envir = .taxodist_cache)
  if (clean) cleaned else lineage
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
      candidates = candidates
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
        )
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
    valid <- vapply(candidate_lineages, function(lineage) {
      !is.null(lineage) && "Biota" %in% lineage
    }, logical(1))
    candidates <- candidates[valid, , drop = FALSE]
    candidate_lineages <- candidate_lineages[valid]

    if (nrow(candidates) == 0L) {
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
        any(grepl(
          paste0("\\b", regexEscape(taxon), "\\b"),
          lineage,
          ignore.case = TRUE
        ))
      }, logical(1))
      if (any(exact)) {
        candidates <- candidates[exact, , drop = FALSE]
        candidate_lineages <- candidate_lineages[exact]
      }
    }

    # Apply the same truncation as get_lineage() so that name-based and
    # resolution-based matrices are identical.
    selected_lineage <- .finalize_lineage(candidate_lineages[[1L]], taxon)
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
      candidates = candidates
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
taxo_from_lineages <- function(lineages, ids = NULL, source = "user-supplied") {
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
  class(result) <- c("taxodist_resolution", "data.frame")
  attr(result, "source") <- source
  attr(result, "source_url") <- NA_character_
  attr(result, "retrieved_at") <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  result
}
