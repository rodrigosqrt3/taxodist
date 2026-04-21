#' @importFrom httr GET content status_code add_headers
#' @importFrom rvest read_html html_nodes html_text html_attr
#' @importFrom stringr str_remove str_remove_all str_extract str_trim
NULL

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
    httr::GET(url, httr::add_headers("User-Agent" = "taxodist R package/0.1")),
    error = function(e) NULL
  )

  if (is.null(res) || httr::status_code(res) != 200) {
    if (verbose) cli::cli_alert_warning("Could not reach Taxonomicon for {taxon}")
    return(NULL)
  }

  page <- rvest::read_html(httr::content(res, "text", encoding = "UTF-8"))
  rows <- rvest::html_nodes(page, "tr")

  for (row in rows) {
    text  <- rvest::html_text(row, trim = TRUE)
    hrefs <- rvest::html_attr(rvest::html_nodes(row, "a"), "href")

    if (grepl("astronomical|planet|Minor planet|comet|asteroid",
              text, ignore.case = TRUE)) next

    tree_links <- hrefs[grepl("TaxonTree\\.aspx.*id=[1-9]", hrefs)]
    if (length(tree_links) > 0) {
      id <- stringr::str_remove(
        stringr::str_extract(tree_links[1], "id=([0-9]+)"), "id="
      )

      # verify this is a biological entry by checking lineage root
      candidate_lin <- get_lineage_by_id(id, clean = TRUE, verbose = FALSE)
      if (is.null(candidate_lin) || !"Biota" %in% candidate_lin) next

      assign(cache_key, id, envir = .taxodist_cache)
      if (verbose) cli::cli_alert_success("Found {taxon} with ID {id}")
      return(id)
    }
  }

  if (verbose) cli::cli_alert_warning("{taxon} not found in Taxonomicon")
  NULL
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
    httr::GET(url, httr::add_headers("User-Agent" = "taxodist R package/0.1")),
    error = function(e) NULL
  )
  if (is.null(res) || httr::status_code(res) != 200) {
    if (verbose) cli::cli_alert_warning("Could not retrieve lineage for ID {taxon_id}")
    return(NULL)
  }
  page  <- rvest::read_html(httr::content(res, "text", encoding = "UTF-8"))
  links <- rvest::html_nodes(page, "a[href*='TaxonTree']")

  # remove navigation/menu links — these have no numeric id parameter
  hrefs_all <- rvest::html_attr(links, "href")
  links <- links[grepl("id=[0-9]", hrefs_all)]
  hrefs_all <- hrefs_all[grepl("id=[0-9]", hrefs_all)]

  # drop genus/species child entries
  texts_raw <- rvest::html_text(links, trim = TRUE)
  is_child <- grepl("^(\\[plesion.*\\] )?Genus |^(\\[plesion.*\\] )?Species |incertae sedis", texts_raw)
  first_child <- which(is_child)[1]
  if (!is.na(first_child)) {
    links <- links[seq_len(first_child - 1)]
  }

  hrefs <- rvest::html_attr(links, "href")
  own_pattern <- paste0("id=", taxon_id, "($|&)")
  own_idx <- which(grepl(own_pattern, hrefs))
  if (length(own_idx) > 0) {
    links <- links[seq_len(max(own_idx))]
  }

  texts <- rvest::html_text(links, trim = TRUE)
  lineage <- texts |>
    stringr::str_remove("^\\[crown\\]\\s+(Clade|Grandorder|Order|Superorder|Infraorder|Suborder|Class|Superclass|Subclass|Infraclass|Family|Superfamily|Subfamily|Tribe|Subtribe|Kingdom|Subkingdom|Infrakingdom|Domain|Superkingdom|Phylum|Subphylum|Genus|Species)?\\s*") |>
    stringr::str_remove("^(Clade |Kingdom |Phylum |Class |Order |Suborder |Infraorder |Parvorder |Grandorder |Magnorder |Cohort |Subcohort |Legion |Family |Subfamily |Tribe |Subtribe |Genus |Species |Subkingdom |Infrakingdom |Superclass |Subclass |Infraclass |Superorder |Superfamily |Domain |Superkingdom )") |>
    stringr::str_remove("\\s+[A-Z][a-z\u00e1\u00e0\u00e2\u00e3\u00e9\u00e8\u00ea\u00ed\u00ef\u00f3\u00f4\u00f5\u00f6\u00fa\u00fc\u00e7].*$") |>
    stringr::str_remove("\\s+[A-Z]\\.[A-Z]\\..*$") |>
    stringr::str_remove("\\s+von.*$") |>
    stringr::str_remove_all("[\u2020\u1D40]") |>
    stringr::str_remove("\\s+\\([A-Z][a-z].*$") |>
    stringr::str_remove("\\s+\\(\\d{4}\\).*$") |>
    stringr::str_remove("\\s+\\[.*$") |>
    stringr::str_remove("\\s+[A-Z]\\.$") |>
    stringr::str_remove("\\s+\\([a-z].*$") |>
    stringr::str_remove("^\".*") |>
    stringr::str_trim()
  bare_ranks <- c(
    "Subphylum", "Infraphylum", "Subfamily", "Suborder",
    "Infraorder", "Superclass", "Subclass", "Superfamily",
    "[crown]", ""
  )
  lineage <- lineage[!lineage %in% bare_ranks]
  lineage <- lineage[lineage != "" & !grepl("^\\s*$", lineage)]
  lineage <- lineage[!grepl("^\"", lineage)]
  lineage <- lineage[!grepl("^Population", lineage)]
  if (clean) {
    phil_nodes <- c("Natura - nature", "actualia - actual entities",
                    "Mundus", "naturalia - natural bodies")
    lineage <- lineage[!lineage %in% phil_nodes]
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
  id <- get_taxonomicon_id(taxon, verbose = verbose)
  if (is.null(id)) return(NULL)
  lineage <- get_lineage_by_id(id, clean = clean, verbose = verbose)
  if (is.null(lineage)) return(NULL)

  if (!grepl("\\s", taxon)) {
    lineage <- lineage[!grepl(" ", lineage)]
    lineage <- lineage[!grepl("^\\[", lineage)]
    # truncate at taxon if present, otherwise append it
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

  if (length(lineage) == 0) return(NULL) # nocov
  lineage
}

