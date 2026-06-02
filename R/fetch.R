#' @importFrom httr GET content status_code add_headers
#' @importFrom rvest read_html html_nodes html_text html_attr
#' @importFrom stringr str_remove str_remove_all str_extract str_trim
NULL

regexEscape <- function(x) gsub("([.\\^$*+?{}\\[\\]|()])", "\\\\\\1", x)

# -- Internal cache -----------------------------------------------------------

.taxodist_cache <- new.env(parent = emptyenv())

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
  rm(list = ls(.taxodist_cache), envir = .taxodist_cache)
  invisible(NULL)
}

#' Save the taxodist lineage cache to disk
#'
#' Serialises the current session cache to an `.rds` file so it can be
#' restored in a future session with [load_cache()]. Useful for
#' reproducibility and for avoiding repeated network requests.
#'
#' @param file Path to the `.rds` file to write.
#'
#' @return Invisibly returns `NULL`.
#' @seealso [load_cache()], [clear_cache()]
#' @export
#' @examples
#' \donttest{
#' tmp <- tempfile(fileext = ".rds")
#' save_cache(tmp)
#' unlink(tmp)
#' }
save_cache <- function(file) {
  data <- as.list(.taxodist_cache)
  saveRDS(data, file = file)
  cli::cli_alert_success("Cache saved to {.file {file}} ({length(data)} entries).")
  invisible(NULL)
}

#' Load a previously saved taxodist cache from disk
#'
#' Restores lineage data saved with [save_cache()] into the current session
#' cache, avoiding network requests for taxa already retrieved in a previous
#' session.
#'
#' @param file Path to an `.rds` file created by [save_cache()].
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
load_cache <- function(file) {
  if (!file.exists(file)) {
    cli::cli_abort("Cache file not found: {.file {file}}")
  }
  data <- readRDS(file)
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
#'   \item{`taxa`}{Character vector of taxa with cached lineages.}
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
  keys <- ls(.taxodist_cache)

  lin_keys <- keys[startsWith(keys, "lin_")]
  id_keys  <- keys[startsWith(keys, "id_")]

  taxa_names <- sub("^lin_", "", lin_keys)

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
      paste0("{.val ", taxa_names, "}"),
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
#'   ranks are also supported.
#' @param verbose Logical. If `TRUE`, prints status messages during retrieval.
#'   Default is `FALSE`.
#'
#' @return A character string containing the Taxonomicon numeric ID, or `NULL`
#'   if the taxon is not found.
#'
#' @details
#' The function queries the static search endpoint at
#' `taxonomicon.taxonomy.nl/TaxonList.aspx` and parses the resulting HTML
#' to extract the taxon ID from the hierarchy link. When multiple matches
#' exist (e.g., a genus name shared with an astronomical object), biological
#' entries are prioritised by filtering for entries annotated as dinosaur,
#' reptile, archosaur, animal, plant, fungus, or bacterium.
#'
#' @seealso [get_lineage()], [taxo_distance()]
#'
#' @export
#' @examples
#' \donttest{
#' get_taxonomicon_id("Tyrannosaurus")   # returns "50841"
#' get_taxonomicon_id("Homo")
#' get_taxonomicon_id("Quercus")
#' }
get_taxonomicon_id <- function(taxon, verbose = FALSE) {
  cache_key <- paste0("id_", taxon)
  if (exists(cache_key, envir = .taxodist_cache)) {
    if (verbose) cli::cli_alert_info("Using cached ID for {taxon}")
    return(get(cache_key, envir = .taxodist_cache))
  }

  if (verbose) cli::cli_alert_info("Searching Taxonomicon for {taxon}...")

  url <- paste0(
    "http://taxonomicon.taxonomy.nl/TaxonList.aspx",
    "?subject=Entity&by=ScientificName&search=",
    utils::URLencode(taxon)
  )

  res <- tryCatch(
    httr::GET(url, httr::add_headers("User-Agent" = "taxodist R package/0.3"), httr::timeout(30)),
    error = function(e) NULL
  )

  if (is.null(res) || httr::status_code(res) != 200) {
    cli::cli_warn(c(
      "!" = "Cannot reach The Taxonomicon server.",
      "x" = "The website (taxonomy.nl) appears to be offline or unreachable.",
      "i" = "Please try again later."
    ))
    return(NULL)
  }

  page <- rvest::read_html(httr::content(res, "text", encoding = "UTF-8"))
  rows <- rvest::html_nodes(page, "tr")
  bio_ids <- list()

  for (row in rows) {
    text <- rvest::html_text(row, trim = TRUE)
    if (grepl("\\bastronomical\\b|\\bplanet\\b|\\bMinor planet\\b|\\bcomet\\b|\\basteroid\\b",
              text, ignore.case = TRUE)) next

    links_nodes <- rvest::html_nodes(row, "a[href*='TaxonTree']")
    if (length(links_nodes) == 0) next

    classes <- rvest::html_attr(links_nodes, "class")
    valid_idx <- which(!is.na(classes) & classes == "Valid")
    if (length(valid_idx) == 0) next

    target_link <- links_nodes[[valid_idx[1]]]

    id <- stringr::str_remove(
      stringr::str_extract(
        rvest::html_attr(target_link, "href"), "id=([0-9]+)"
      ), "id="
    )
    if (is.na(id)) next

    text_entry <- trimws(gsub("\\s+", " ", text))
    text_entry <- stringr::str_remove(text_entry, "^N\\|T\\|P\\|R\\|B\\|L\\s*")

    candidate_lin <- get_lineage_by_id(id, clean = TRUE, verbose = FALSE)
    if (is.null(candidate_lin) || !"Biota" %in% candidate_lin) next

    bio_ids <- c(bio_ids, list(list(id = id, text = text_entry)))
  }

  if (length(bio_ids) > 1) {
    matched <- Filter(function(x) {
      lin <- get_lineage_by_id(x$id, clean = TRUE, verbose = FALSE)
      !is.null(lin) && any(grepl(paste0("\\b", regexEscape(taxon), "\\b"), lin, ignore.case = TRUE))
    }, bio_ids)
    # nocov start
    if (length(matched) > 0) bio_ids <- matched
    # nocov end
  }

  if (length(bio_ids) == 0L) {
    if (verbose) cli::cli_alert_warning("{taxon} not found in Taxonomicon")
    return(NULL)
  }

  unique_ids <- unique(sapply(bio_ids, function(x) x$id))
  bio_ids <- lapply(unique_ids, function(uid) {
    matches <- bio_ids[sapply(bio_ids, function(x) x$id == uid)]
    matches[[1]]
  })

  if (length(bio_ids) > 1L) {
    warn_msg <- c(
      "!" = "Multiple valid biological entries found for {.val {taxon}}.",
      "i" = paste0("Using: ", bio_ids[[1]]$text, " (ID: ", bio_ids[[1]]$id, ")"),
      "i" = paste0("To use a different entry, pass its numeric ID directly, e.g. `get_lineage(\"", bio_ids[[2]]$id, "\")`."),
      "i" = "Other available IDs:"
    )
    for (i in 2:length(bio_ids)) {
      warn_msg <- c(warn_msg, "*" = paste0("ID ", bio_ids[[i]]$id, ": ", bio_ids[[i]]$text))
    }
    cli::cli_warn(warn_msg)
  }

  id <- bio_ids[[1]]$id
  assign(cache_key, id, envir = .taxodist_cache)
  if (verbose) cli::cli_alert_success("Found {taxon} with ID {id}")
  return(id)
}

# -- Lineage retrieval --------------------------------------------------------

#' Retrieve the full taxonomic lineage of a taxon
#'
#' Given a Taxonomicon numeric ID, retrieves and parses the complete
#' hierarchical lineage from root (Natura) to the taxon itself. The lineage
#' is returned as a character vector ordered from root to tip.
#'
#' @param taxon_id A numeric or character string giving the Taxonomicon ID.
#'   Obtain this with [get_taxonomicon_id()].
#' @param clean Logical. If `TRUE` (default), removes philosophical root nodes
#'   above `Biota` (i.e., Natura, actualia, Mundus, naturalia) and strips
#'   dagger and superscript markers from names.
#' @param verbose Logical. If `TRUE`, prints status messages. Default `FALSE`.
#'
#' @return A character vector of clade names from root to tip, or `NULL` if
#'   retrieval fails.
#'
#' @details
#' Lineage data is sourced from The Taxonomicon, which is based on
#' Systema Naturae 2000 (Brands, S.J., 1989 onwards). The depth of lineages
#' in The Taxonomicon substantially exceeds that of other programmatic sources
#' such as the Open Tree of Life, particularly for well-studied clades such
#' as Dinosauria, where intermediate clades at the level of superfamilies,
#' tribes, and named subclades are fully resolved.
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
  if (is.null(taxon_id) || is.na(taxon_id) || taxon_id == "" || !grepl("^[0-9]+$", as.character(taxon_id))) {
    return(NULL)
  }
  cache_key <- paste0("lin_", taxon_id)
  if (exists(cache_key, envir = .taxodist_cache)) {
    if (verbose) cli::cli_alert_info("Using cached lineage for ID {taxon_id}")
    return(get(cache_key, envir = .taxodist_cache))
  }
  url <- paste0(
    "http://taxonomicon.taxonomy.nl/TaxonTree.aspx?id=",
    taxon_id, "&src=0"
  )
  res <- tryCatch(
    httr::GET(
      url,
      httr::add_headers("User-Agent" = "taxodist R package/0.3"),
      httr::timeout(30)
    ),
    error = function(e) NULL
  )
  if (is.null(res) || httr::status_code(res) != 200) {
    if (verbose) cli::cli_alert_warning("Could not retrieve lineage for ID {taxon_id}")
    return(NULL)
  }
  page  <- rvest::read_html(httr::content(res, "text", encoding = "UTF-8"))
  regexEscape <- function(x) gsub("([.\\^$*+?{}\\[\\]|()])", "\\\\\\1", x)

  current_name <- rvest::html_text(
    rvest::html_node(page, "#ctl00_divSubject b"), trim = TRUE
  )
  content_node <- rvest::html_node(page, "#divPageContent")

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
    links <- rvest::html_nodes(page, "a[href*='TaxonTree']")
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
  if (clean) {
    idx_biota <- which(lineage == "Biota")
    if (length(idx_biota) > 0) {
      lineage <- lineage[idx_biota[1]:length(lineage)]
    }
  }
  if (length(lineage) == 0) return(NULL)
  assign(cache_key, lineage, envir = .taxodist_cache)
  lineage
}

#' Retrieve the full taxonomic lineage of a taxon by name
#'
#' A convenience wrapper that combines [get_taxonomicon_id()] and
#' [get_lineage_by_id()] into a single call. Given a taxon name, returns
#' its complete lineage from root to tip.
#'
#' @param taxon A character string giving the taxon name.
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
  is_id <- grepl("^[0-9]+$", as.character(taxon))

  if (is_id) {
    id <- as.character(taxon)
  } else {
    id <- get_taxonomicon_id(taxon, verbose = verbose)
  }

  if (is.null(id)) return(NULL)

  lineage <- get_lineage_by_id(id, clean = clean, verbose = verbose)
  if (is.null(lineage)) return(NULL)

  if (!is_id) {
    if (!grepl("\\s", taxon)) {
      lineage <- lineage[!grepl(" ", lineage)]
      lineage <- lineage[!grepl("^\\[", lineage)]
      target_idx <- which(lineage == taxon)
      if (length(target_idx) > 0) {
        lineage <- lineage[seq_len(target_idx[length(target_idx)])]
      } else {
        lineage <- c(lineage, taxon)
      }
    } else {
      lineage <- lineage[!grepl(" ", lineage) | lineage == taxon]
      target_idx <- which(lineage == taxon)
      if (length(target_idx) > 0) {
        lineage <- lineage[seq_len(target_idx[1])]
      } else {
        lineage <- c(lineage, taxon) # nocov
      }
    }
  }

  if (length(lineage) == 0) return(NULL) # nocov
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
#' @seealso [get_lineage()], [taxo_distance()]
#' @export
#' @examples
#' \donttest{
#' taxo_search("Bacteria")
#' taxo_search("Nereis")
#' taxo_search("Tyrannosaurus")
#' }
taxo_search <- function(taxon, verbose = FALSE) {
  if (verbose) cli::cli_alert_info("Searching Taxonomicon for {.val {taxon}}...")

  url <- paste0(
    "http://taxonomicon.taxonomy.nl/TaxonList.aspx",
    "?subject=Entity&by=ScientificName&search=",
    utils::URLencode(taxon)
  )

  res <- tryCatch(
    httr::GET(url, httr::add_headers("User-Agent" = "taxodist R package/0.3"), httr::timeout(30)),
    error = function(e) NULL
  )

  if (is.null(res) || httr::status_code(res) != 200) {
    if (verbose) cli::cli_alert_warning("Could not reach Taxonomicon")
    return(NULL)
  }

  page <- rvest::read_html(httr::content(res, "text", encoding = "UTF-8"))
  rows <- rvest::html_nodes(page, "tr")

  results <- list()
  for (row in rows) {
    text <- rvest::html_text(row, trim = TRUE)
    if (grepl("\\bastronomical\\b|\\bplanet\\b|\\bMinor planet\\b|\\bcomet\\b|\\basteroid\\b", text, ignore.case = TRUE)) next

    links <- rvest::html_nodes(row, "a[href*='TaxonTree']")
    if (length(links) == 0) next

    classes <- rvest::html_attr(links, "class")
    valid_idx <- which(!is.na(classes) & classes == "Valid")
    if (length(valid_idx) == 0) next

    id <- stringr::str_remove(stringr::str_extract(rvest::html_attr(links[[valid_idx[1]]], "href"), "id=([0-9]+)"), "id=")
    if (is.na(id)) next

    text_entry <- trimws(gsub("\\s+", " ", text))
    text_entry <- stringr::str_remove(text_entry, "^N\\|T\\|P\\|R\\|B\\|L\\s*")

    results <- c(results, list(data.frame(id = id, name = text_entry, stringsAsFactors = FALSE)))
  }

  if (length(results) == 0) {
    if (verbose) cli::cli_alert_warning("No matches found.")
    return(NULL)
  }

  df <- do.call(rbind, results)
  df <- df[!duplicated(df$id), ]
  rownames(df) <- NULL

  if (verbose) cli::cli_alert_success("Found {nrow(df)} entries.")
  return(df)
}

