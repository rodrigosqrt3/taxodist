# Taxonomic sources other than The Taxonomicon.
#
# Every source resolves names to records with the same structure as
# taxo_resolve() uses for The Taxonomicon, so the resulting
# `taxodist_resolution` works with every matrix, tree, bundle, and
# distinctness function. Identifiers carry a source prefix ("ott:770309",
# "ncbi:9606") so they are never confused with Taxonomicon IDs.

.source_ids <- c("taxonomicon", "ott", "ncbi")

.ott_api <- "https://api.opentreeoflife.org/v3"
.ncbi_eutils <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

#' List the taxonomic sources supported by taxodist
#'
#' Describes the classifications from which [taxo_resolve()] can retrieve
#' lineages, with their licence and preferred citation. Distances are only
#' comparable within one source; use [compare_sources()] to compare sources.
#'
#' @return A data frame with one row per source and columns `id` (the value
#'   accepted by the `source` argument of [taxo_resolve()]), `name`, `url`,
#'   `id_prefix` (prefix of the identifiers stored in resolutions), `ranks`
#'   (whether unranked clades are included), `licence`, and `citation`.
#'
#' @details
#' \describe{
#'   \item{`"taxonomicon"`}{The Taxonomicon (Brands 1989--present), a
#'     single curated classification with very deep lineages, many unranked
#'     clades, and extensive fossil coverage. Retrieved from its web pages.}
#'   \item{`"ott"`}{The Open Tree of Life reference taxonomy (Rees and
#'     Cranston 2017), a synthesis of NCBI, GBIF, WoRMS, IRMNG, and other
#'     taxonomies, retrieved from the Open Tree of Life API. Its version is
#'     recorded in the resolution.}
#'   \item{`"ncbi"`}{The NCBI Taxonomy database (Schoch et al. 2020),
#'     covering organisms with sequence data, retrieved from the NCBI
#'     E-utilities. NCBI allows three requests per second, or ten with an
#'     API key set in the `NCBI_API_KEY` environment variable or the
#'     `taxodist.ncbi_api_key` option; requests are spaced by at least 0.4 s
#'     (0.12 s with a key), and refused requests (HTTP 429) are retried.}
#' }
#' Lineages are ordered from the root of each source (`"Biota"`, `"life"`,
#' and `"cellular organisms"` respectively) to the taxon. Clade names differ
#' between sources: for example, animals are `"Animalia"` in The Taxonomicon
#' but `"Metazoa"` in the other two, which matters for the `within` argument
#' of [taxo_resolve()].
#'
#' @references
#' Brands, S. J. (ed.) 1989--present. The Taxonomicon. Universal Taxonomic
#' Services, Zwaag, The Netherlands. <http://taxonomicon.taxonomy.nl>
#'
#' Rees, J. A. and Cranston, K. (2017). Automated assembly of a reference
#' taxonomy for phylogenetic data synthesis. *Biodiversity Data Journal*,
#' 5, e12581. \doi{10.3897/BDJ.5.e12581}
#'
#' Schoch, C. L. et al. (2020). NCBI Taxonomy: a comprehensive update on
#' curation, resources and tools. *Database*, 2020, baaa062.
#' \doi{10.1093/database/baaa062}
#'
#' @seealso [taxo_resolve()], [compare_sources()]
#' @export
#' @examples
#' taxo_sources()
taxo_sources <- function() {
  data.frame(
    id = .source_ids,
    name = c("The Taxonomicon", "Open Tree Taxonomy", "NCBI Taxonomy"),
    url = c(
      "http://taxonomicon.taxonomy.nl",
      "https://tree.opentreeoflife.org/about/taxonomy-version",
      "https://www.ncbi.nlm.nih.gov/taxonomy"
    ),
    id_prefix = c("", "ott:", "ncbi:"),
    ranks = c(
      "ranked and unranked clades",
      "ranked and unranked clades",
      "ranked and unranked clades"
    ),
    licence = c(
      "All rights reserved; cite the source",
      "CC0 (Open Tree of Life); underlying sources keep their terms",
      "Public domain (US government work)"
    ),
    citation = c(
      "Brands, S.J. (ed.) 1989-present. The Taxonomicon. Universal Taxonomic Services, Zwaag, The Netherlands.",
      "Rees, J.A. & Cranston, K. (2017). Automated assembly of a reference taxonomy for phylogenetic data synthesis. Biodiversity Data Journal 5: e12581.",
      "Schoch, C.L. et al. (2020). NCBI Taxonomy: a comprehensive update on curation, resources and tools. Database 2020: baaa062."
    ),
    stringsAsFactors = FALSE
  )
}

.match_source <- function(source) {
  if (!is.character(source) || length(source) != 1L || is.na(source)) {
    cli::cli_abort(
      "{.arg source} must be a single string, one of {.or {.val {(.source_ids)}}}."
    )
  }
  if (!source %in% .source_ids) {
    cli::cli_abort(c(
      "x" = "Unknown taxonomic source {.val {source}}.",
      "i" = "Use one of {.or {.val {(.source_ids)}}}; see {.fn taxo_sources}."
    ))
  }
  source
}

# -- HTTP ------------------------------------------------------------------

# Minimum interval between requests to each source. NCBI allows at most
# three requests per second without an API key and ten with one.
.source_interval <- function(source) {
  interval <- getOption("taxodist.min_interval", 0.25)
  if (!is.numeric(interval) || length(interval) != 1L || is.na(interval) ||
      interval < 0) {
    interval <- 0
  }
  if (identical(source, "ncbi")) {
    # A margin over 1/3 s (1/10 s with a key): NCBI counts requests over a
    # sliding window, and bursts at exactly the limit are refused with 429.
    floor <- if (nzchar(.ncbi_api_key())) 0.12 else 0.4
    interval <- max(interval, floor)
  }
  interval
}

.ncbi_api_key <- function() {
  key <- getOption("taxodist.ncbi_api_key", Sys.getenv("NCBI_API_KEY", ""))
  if (!is.character(key) || length(key) != 1L || is.na(key)) "" else key
}

# Performs one request and returns list(status, text). `status` is the HTTP
# status code, or NA when the request failed without a response.
.source_request <- function(source, method, url, body = NULL) {
  .throttle(source, .source_interval(source))
  args <- list(
    method, url,
    httr::add_headers("User-Agent" = .taxodist_user_agent()),
    httr::timeout(30),
    # Client errors are final, except 429 (too many requests), which is
    # retried after the backoff pause.
    times = .max_tries(), pause_base = 1, pause_cap = 5,
    terminate_on = setdiff(400:499, 429L), quiet = TRUE
  )
  if (!is.null(body)) {
    args$body <- jsonlite::toJSON(body, auto_unbox = TRUE)
    args <- c(args, list(httr::content_type_json()))
  }
  res <- tryCatch(do.call(httr::RETRY, args), error = function(e) NULL)
  if (is.null(res)) return(list(status = NA_integer_, text = NULL))
  list(
    status = as.integer(httr::status_code(res)),
    text = httr::content(res, "text", encoding = "UTF-8")
  )
}

.parse_json <- function(text) {
  if (is.null(text)) return(NULL)
  tryCatch(jsonlite::fromJSON(text, simplifyVector = FALSE), error = function(e) NULL)
}

# Rank labels as used by The Taxonomicon: capitalised, with unranked nodes
# (other than explicit clades) as NA.
.normalize_source_rank <- function(rank) {
  rank <- trimws(tolower(as.character(rank)))
  unranked <- is.na(rank) | !nzchar(rank) |
    rank %in% c("no rank", "no rank - terminal", "cellular root", "acellular root")
  out <- paste0(toupper(substr(rank, 1L, 1L)), substring(rank, 2L))
  out[unranked] <- NA_character_
  out
}

.source_cache_get <- function(key) {
  if (exists(key, envir = .taxodist_cache, inherits = FALSE)) {
    get(key, envir = .taxodist_cache, inherits = FALSE)
  } else {
    NULL
  }
}

.source_cache_set <- function(key, value) {
  assign(key, value, envir = .taxodist_cache)
  invisible(value)
}

# -- Open Tree Taxonomy ----------------------------------------------------

# Version label of the Open Tree Taxonomy served by the API, e.g.
# "Open Tree Taxonomy 3.7 (ott3.7draft3)".
.ott_about <- function() {
  cached <- .source_cache_get("src_ott_about")
  if (!is.null(cached)) return(cached)
  res <- .source_request("ott", "POST", paste0(.ott_api, "/taxonomy/about"), body = list())
  data <- if (identical(res$status, 200L)) .parse_json(res$text) else NULL
  if (is.null(data)) {
    return(list(name = "Open Tree Taxonomy",
                url = "https://tree.opentreeoflife.org/about/taxonomy-version"))
  }
  about <- list(
    name = paste0("Open Tree Taxonomy ", data$version, " (", data$source, ")"),
    url = if (is.null(data$weburl)) {
      "https://tree.opentreeoflife.org/about/taxonomy-version"
    } else {
      data$weburl
    }
  )
  .source_cache_set("src_ott_about", about)
}

# Candidates for each name: a data frame with columns id, name, and synonym,
# an empty data frame when the name is unknown, or NULL when the request
# failed. Names are sent in batches.
.ott_match_names <- function(names, verbose = FALSE) {
  out <- stats::setNames(vector("list", length(names)), names)
  todo <- character(0)
  for (name in names) {
    cached <- .source_cache_get(paste0("src_ott_match_", tolower(name)))
    if (is.null(cached)) todo <- c(todo, name) else out[[name]] <- cached
  }
  batches <- split(todo, ceiling(seq_along(todo) / 250))
  for (batch in batches) {
    if (verbose) cli::cli_alert_info("Matching {length(batch)} name{?s} in the Open Tree Taxonomy...")
    res <- .source_request(
      "ott", "POST", paste0(.ott_api, "/tnrs/match_names"),
      body = list(names = I(batch), do_approximate_matching = FALSE)
    )
    data <- if (identical(res$status, 200L)) .parse_json(res$text) else NULL
    if (is.null(data)) next
    by_query <- list()
    for (result in data$results) by_query[[tolower(result$name)]] <- result$matches
    for (name in batch) {
      matches <- by_query[[tolower(name)]]
      candidates <- data.frame(
        id = vapply(matches, function(m) paste0("ott:", m$taxon$ott_id), character(1)),
        name = vapply(matches, function(m) m$taxon$name, character(1)),
        synonym = vapply(matches, function(m) isTRUE(m$is_synonym), logical(1)),
        stringsAsFactors = FALSE
      )
      .source_cache_set(paste0("src_ott_match_", tolower(name)), candidates)
      out[[name]] <- candidates
    }
  }
  out
}

# Lineage of an OTT taxon from the root ("life") to the taxon, with ranks
# as attribute; NULL when it cannot be retrieved.
.ott_lineage <- function(id, verbose = FALSE) {
  number <- sub("^ott:", "", id)
  key <- paste0("src_ott_lin_", number)
  cached <- .source_cache_get(key)
  if (!is.null(cached)) return(cached)
  if (verbose) cli::cli_alert_info("Retrieving lineage of ott:{number}...")
  res <- .source_request(
    "ott", "POST", paste0(.ott_api, "/taxonomy/taxon_info"),
    body = list(ott_id = as.numeric(number), include_lineage = TRUE)
  )
  data <- if (identical(res$status, 200L)) .parse_json(res$text) else NULL
  if (is.null(data) || is.null(data$name)) return(NULL)
  ancestors <- rev(data$lineage)
  lineage <- c(
    vapply(ancestors, function(node) node$name, character(1)),
    data$name
  )
  ranks <- .normalize_source_rank(c(
    vapply(ancestors, function(node) {
      if (is.null(node$rank)) NA_character_ else node$rank
    }, character(1)),
    if (is.null(data$rank)) NA_character_ else data$rank
  ))
  lineage <- structure(lineage, ranks = ranks)
  .source_cache_set(key, lineage)
}

# -- NCBI Taxonomy ---------------------------------------------------------

.ncbi_url <- function(endpoint, query) {
  key <- .ncbi_api_key()
  query <- c(query, tool = "taxodist", if (nzchar(key)) c(api_key = key))
  paste0(
    .ncbi_eutils, "/", endpoint, "?",
    paste(names(query), utils::URLencode(query, reserved = TRUE), sep = "=", collapse = "&")
  )
}

# NCBI identifiers matching a name: by scientific name, then by any name
# (synonyms, common names). character(0) when unknown; NULL on failure.
.ncbi_search <- function(name, verbose = FALSE) {
  key <- paste0("src_ncbi_search_", tolower(name))
  cached <- .source_cache_get(key)
  if (!is.null(cached)) return(cached)
  if (verbose) cli::cli_alert_info("Searching NCBI Taxonomy for {name}...")
  for (field in c("Scientific Name", "All Names")) {
    term <- paste0("\"", gsub("\"", "", name), "\"[", field, "]")
    res <- .source_request(
      "ncbi", "GET",
      .ncbi_url("esearch.fcgi", c(db = "taxonomy", term = term, retmode = "json", retmax = "50"))
    )
    data <- if (identical(res$status, 200L)) .parse_json(res$text) else NULL
    if (is.null(data) || is.null(data$esearchresult)) return(NULL)
    ids <- unlist(data$esearchresult$idlist, use.names = FALSE)
    if (length(ids) > 0L) {
      ids <- paste0("ncbi:", ids)
      return(.source_cache_set(key, ids))
    }
  }
  .source_cache_set(key, character(0))
}

# Retrieves and caches the lineages of NCBI taxa, in batches.
.ncbi_fetch <- function(ids, verbose = FALSE) {
  numbers <- unique(sub("^ncbi:", "", ids))
  todo <- numbers[!vapply(numbers, function(n) {
    exists(paste0("src_ncbi_lin_", n), envir = .taxodist_cache, inherits = FALSE)
  }, logical(1))]
  for (batch in split(todo, ceiling(seq_along(todo) / 200))) {
    if (verbose) cli::cli_alert_info("Retrieving {length(batch)} lineage{?s} from NCBI Taxonomy...")
    res <- .source_request(
      "ncbi", "GET",
      .ncbi_url("efetch.fcgi", c(db = "taxonomy", id = paste(batch, collapse = ","), retmode = "xml"))
    )
    if (!identical(res$status, 200L)) next
    .ncbi_parse_taxa(res$text)
  }
  invisible(NULL)
}

# Parses an efetch TaxaSet and caches the lineage of every taxon in it.
.ncbi_parse_taxa <- function(text) {
  doc <- tryCatch(xml2::read_xml(text), error = function(e) NULL)
  if (is.null(doc)) return(invisible(NULL))
  for (taxon in xml2::xml_find_all(doc, "/TaxaSet/Taxon")) {
    id <- xml2::xml_text(xml2::xml_find_first(taxon, "./TaxId"))
    name <- xml2::xml_text(xml2::xml_find_first(taxon, "./ScientificName"))
    rank <- xml2::xml_text(xml2::xml_find_first(taxon, "./Rank"))
    ancestors <- xml2::xml_find_all(taxon, "./LineageEx/Taxon")
    lineage <- c(
      xml2::xml_text(xml2::xml_find_first(ancestors, "./ScientificName")),
      name
    )
    ranks <- .normalize_source_rank(c(
      xml2::xml_text(xml2::xml_find_first(ancestors, "./Rank")),
      rank
    ))
    .source_cache_set(paste0("src_ncbi_lin_", id), structure(lineage, ranks = ranks))
    # A merged identifier answers with the taxon it was merged into.
    for (old in xml2::xml_text(xml2::xml_find_all(taxon, "./AkaTaxIds/TaxId"))) {
      .source_cache_set(paste0("src_ncbi_lin_", old), structure(lineage, ranks = ranks))
    }
  }
  invisible(NULL)
}

.ncbi_lineage <- function(id, verbose = FALSE) {
  number <- sub("^ncbi:", "", id)
  key <- paste0("src_ncbi_lin_", number)
  if (!exists(key, envir = .taxodist_cache, inherits = FALSE)) {
    .ncbi_fetch(number, verbose = verbose)
  }
  .source_cache_get(key)
}

# -- Resolution ------------------------------------------------------------

.source_lineage <- function(source, id, verbose = FALSE) {
  switch(source,
    ott = .ott_lineage(id, verbose = verbose),
    ncbi = .ncbi_lineage(id, verbose = verbose)
  )
}

.source_id_pattern <- function(source) {
  paste0("^(", source, ":)?[0-9]+$")
}

# Resolves unique, normalized names in an external source. Returns one
# record per name, with the fields used by taxo_resolve().
.resolve_external <- function(names, source, within_clade, verbose, progress,
                              failed) {
  is_id <- grepl(.source_id_pattern(source), names)
  search_names <- names[!is_id]

  candidates_by_name <- switch(source,
    ott = .ott_match_names(search_names, verbose = verbose),
    ncbi = {
      found <- lapply(search_names, .ncbi_search, verbose = verbose)
      ids <- unlist(found, use.names = FALSE)
      if (length(ids)) .ncbi_fetch(ids, verbose = verbose)
      stats::setNames(lapply(found, function(ids) {
        if (is.null(ids)) return(NULL)
        lineages <- lapply(ids, function(id) .source_cache_get(paste0("src_ncbi_lin_", sub("^ncbi:", "", id))))
        data.frame(
          id = ids,
          name = vapply(lineages, function(l) if (is.null(l)) NA_character_ else utils::tail(l, 1L), character(1)),
          synonym = rep(FALSE, length(ids)),
          stringsAsFactors = FALSE
        )
      }), search_names)
    }
  )

  progress_id <- NULL
  if (progress && length(names) > 0L) {
    progress_id <- cli::cli_progress_bar("Resolving taxa", total = length(names))
  }
  records <- lapply(seq_along(names), function(k) {
    name <- names[[k]]
    record <- if (is_id[[k]]) {
      .resolve_external_id(name, source, verbose, failed)
    } else {
      .resolve_external_name(name, candidates_by_name[[name]], source,
                             within_clade, verbose, failed)
    }
    if (progress) cli::cli_progress_update(id = progress_id)
    record
  })
  if (progress && length(names) > 0L) cli::cli_progress_done(id = progress_id)
  records
}

.resolve_external_id <- function(input, source, verbose, failed) {
  id <- if (startsWith(input, paste0(source, ":"))) input else paste0(source, ":", input)
  lineage <- .source_lineage(source, id, verbose = verbose)
  if (is.null(lineage)) return(failed(input, "retrieval_error"))
  name <- utils::tail(lineage, 1L)
  list(
    input = input,
    resolved_name = name,
    id = id,
    status = "resolved",
    n_candidates = 1L,
    lineage_depth = length(lineage),
    lineage = as.character(lineage),
    candidates = data.frame(id = id, name = name, stringsAsFactors = FALSE),
    ranks = attr(lineage, "ranks", exact = TRUE)
  )
}

.resolve_external_name <- function(input, candidates, source, within_clade,
                                   verbose, failed) {
  if (is.null(candidates)) return(failed(input, "retrieval_error"))
  # NCBI candidates have no name when their lineage could not be retrieved.
  named <- !is.na(candidates$name)
  if (nrow(candidates) > 0L && !any(named)) return(failed(input, "retrieval_error"))
  candidates <- candidates[named, , drop = FALSE]
  if (nrow(candidates) == 0L) return(failed(input, "unresolved"))

  # Exact matches of the queried name take precedence over matches through
  # synonyms or other names, as the lineage match does for The Taxonomicon.
  exact <- tolower(candidates$name) == tolower(input)
  if (any(exact)) candidates <- candidates[exact, , drop = FALSE]
  # NCBI returns identifiers in no particular order; lower (older)
  # identifiers come first, which is deterministic and usually puts the
  # best-known taxon first.
  if (identical(source, "ncbi")) {
    candidates <- candidates[order(as.numeric(sub("^ncbi:", "", candidates$id))), , drop = FALSE]
  }
  candidates <- candidates[order(candidates$synonym), , drop = FALSE]

  needs_all <- !is.null(within_clade) || nrow(candidates) > 1L
  lineages <- vector("list", nrow(candidates))
  for (i in seq_len(nrow(candidates))) {
    if (needs_all || i == 1L) {
      # `[<-` with list() keeps NULL (a failed retrieval) as an element.
      lineages[i] <- list(.source_lineage(source, candidates$id[[i]], verbose = verbose))
    }
  }
  if (nrow(candidates) > 1L) {
    keep <- .outermost_candidates(lineages)
    candidates <- candidates[keep, , drop = FALSE]
    lineages <- lineages[keep]
  }
  if (!is.null(within_clade)) {
    retrieved <- !vapply(lineages, is.null, logical(1))
    if (!any(retrieved)) {
      return(failed(input, "retrieval_error", candidates[, c("id", "name")]))
    }
    inside <- retrieved & vapply(lineages, function(lineage) {
      !is.null(lineage) && within_clade %in% tolower(lineage)
    }, logical(1))
    if (!any(inside)) return(failed(input, "unresolved"))
    candidates <- candidates[inside, , drop = FALSE]
    lineages <- lineages[inside]
  }
  lineage <- lineages[[1L]]
  if (is.null(lineage)) {
    return(failed(input, "retrieval_error", candidates[, c("id", "name")]))
  }
  n <- nrow(candidates)
  candidates <- candidates[, c("id", "name"), drop = FALSE]
  rownames(candidates) <- NULL
  list(
    input = input,
    resolved_name = utils::tail(lineage, 1L),
    id = candidates$id[[1L]],
    status = if (n > 1L) "ambiguous" else "resolved",
    n_candidates = n,
    lineage_depth = length(lineage),
    lineage = as.character(lineage),
    candidates = candidates,
    ranks = attr(lineage, "ranks", exact = TRUE)
  )
}

# Candidates nested in another candidate of the same name, such as the
# nominotypical subgenus Mus (Mus) inside the genus Mus, are not homonyms:
# only the outermost is kept. A candidate is dropped when the lineage of
# another candidate is a proper prefix of its own. Candidates without a
# lineage are kept.
.outermost_candidates <- function(lineages) {
  n <- length(lineages)
  vapply(seq_len(n), function(i) {
    li <- lineages[[i]]
    if (is.null(li)) return(TRUE)
    !any(vapply(seq_len(n), function(j) {
      lj <- lineages[[j]]
      j != i && !is.null(lj) && length(lj) < length(li) &&
        identical(as.character(li[seq_along(lj)]), as.character(lj))
    }, logical(1)))
  }, logical(1))
}

.source_metadata <- function(source) {
  switch(source,
    taxonomicon = list(name = "The Taxonomicon", url = "http://taxonomicon.taxonomy.nl"),
    ott = .ott_about(),
    ncbi = list(name = "NCBI Taxonomy", url = "https://www.ncbi.nlm.nih.gov/taxonomy")
  )
}

# -- Comparison of sources -------------------------------------------------

#' Compare taxonomic distances across sources
#'
#' Resolves the same taxa in several taxonomic sources and compares the
#' resulting classifications: the depth of every lineage in each source and
#' the agreement between the distance matrices, measured by rank correlation
#' and a Mantel permutation test. This shows how much an analysis depends on
#' the classification it is based on.
#'
#' @param taxa A character vector of taxon names (identifiers are specific to
#'   one source and are not accepted).
#' @param sources Character vector of at least two source identifiers; see
#'   [taxo_sources()].
#' @param metric Distance metric; see [taxo_distance()].
#' @param ambiguity How ambiguous names are handled in each source; see
#'   [taxo_resolve()]. Defaults to `"first"`, so ambiguity is reported in the
#'   result rather than as warnings.
#' @param within Optional clade restricting the candidates, either one name
#'   used for every source or a named character vector with one name per
#'   source (clade names differ between sources, e.g. `"Animalia"` and
#'   `"Metazoa"`).
#' @param permutations Number of permutations of the Mantel test; `0` skips
#'   the test.
#' @param verbose,progress Passed to [taxo_resolve()].
#'
#' @return An object of class `taxodist_source_comparison`, a list with:
#' \describe{
#'   \item{`resolutions`}{Named list of `taxodist_resolution` objects.}
#'   \item{`matrices`}{Named list of distance matrices restricted to the taxa
#'     resolved in every source (`shared`).}
#'   \item{`shared`}{The taxa resolved in every source.}
#'   \item{`taxa`}{Data frame with one row per input and, for every source,
#'     its status, resolved identifier, and lineage depth.}
#'   \item{`agreement`}{Data frame with one row per pair of sources: the
#'     number of taxa, the Spearman correlation between their distances, and
#'     the Mantel permutation p-value (`NA` when not computed).}
#'   \item{`metric`, `permutations`}{The settings used.}
#' }
#'
#' @details
#' The Mantel test permutes the taxa of one matrix and reports the
#' proportion of permutations whose correlation is at least the observed
#' one, counting the observed arrangement, i.e. `(1 + count) / (1 +
#' permutations)`. Lineage depths are not comparable across sources in
#' absolute terms, since each classification resolves different parts of the
#' tree of life to different degrees; correlations of distances are.
#'
#' @seealso [taxo_sources()], [taxo_resolve()]
#' @export
#' @examplesIf taxonomicon_available()
#' \donttest{
#' cmp <- compare_sources(
#'   c("Homo", "Pan", "Mus", "Gallus", "Danio", "Drosophila"),
#'   sources = c("taxonomicon", "ott", "ncbi")
#' )
#' cmp
#' }
compare_sources <- function(taxa,
                            sources = c("taxonomicon", "ott", "ncbi"),
                            metric = "reciprocal",
                            ambiguity = "first",
                            within = NULL,
                            permutations = 999L,
                            verbose = FALSE,
                            progress = TRUE) {
  metric <- .match_metric(metric)
  if (!is.character(sources) || length(sources) < 2L || anyNA(sources) ||
      anyDuplicated(sources)) {
    cli::cli_abort("{.arg sources} must name at least two different sources.")
  }
  for (source in sources) .match_source(source)
  if (!is.character(taxa) || length(taxa) < 3L) {
    cli::cli_abort("{.arg taxa} must be a character vector of at least three names.")
  }
  if (!is.numeric(permutations) || length(permutations) != 1L ||
      is.na(permutations) || permutations < 0) {
    cli::cli_abort("{.arg permutations} must be a single non-negative number.")
  }
  permutations <- as.integer(permutations)
  within_for <- function(source) {
    if (is.null(within)) return(NULL)
    if (!is.null(names(within))) {
      if (!source %in% names(within)) return(NULL)
      return(unname(within[[source]]))
    }
    if (length(within) != 1L) {
      cli::cli_abort("{.arg within} must be one clade name or a vector named by source.")
    }
    within
  }

  resolutions <- stats::setNames(lapply(sources, function(source) {
    taxo_resolve(
      taxa, ambiguity = ambiguity, verbose = verbose, progress = progress,
      within = within_for(source), source = source
    )
  }), sources)

  usable <- function(res) res$status %in% c("resolved", "ambiguous")
  shared_rows <- Reduce(`&`, lapply(resolutions, usable))
  # Identical inputs are compared once.
  keep <- shared_rows & !duplicated(taxa)
  shared <- taxa[keep]

  matrices <- lapply(resolutions, function(res) {
    sub <- res[keep, , drop = FALSE]
    distance_matrix(sub, progress = FALSE, metric = metric)
  })

  taxa_table <- data.frame(input = taxa, stringsAsFactors = FALSE)
  for (source in sources) {
    res <- resolutions[[source]]
    taxa_table[[paste0("status_", source)]] <- res$status
    taxa_table[[paste0("id_", source)]] <- res$id
    taxa_table[[paste0("depth_", source)]] <- res$lineage_depth
  }

  pairs <- utils::combn(sources, 2L, simplify = FALSE)
  agreement <- do.call(rbind, lapply(pairs, function(pair) {
    a <- matrices[[pair[[1L]]]]
    b <- matrices[[pair[[2L]]]]
    stats <- .mantel(a, b, permutations)
    data.frame(
      source_a = pair[[1L]],
      source_b = pair[[2L]],
      n_taxa = length(shared),
      spearman = stats$r,
      p_value = stats$p,
      stringsAsFactors = FALSE
    )
  }))

  structure(
    list(
      resolutions = resolutions,
      matrices = matrices,
      shared = shared,
      taxa = taxa_table,
      agreement = agreement,
      metric = metric,
      permutations = permutations
    ),
    class = "taxodist_source_comparison"
  )
}

# Spearman Mantel statistic between two `dist` objects over the same taxa,
# with a one-sided permutation p-value. Pairs with a non-finite distance in
# either matrix are excluded.
.mantel <- function(a, b, permutations) {
  x <- as.vector(a)
  y <- as.vector(b)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L) return(list(r = NA_real_, p = NA_real_))
  r <- suppressWarnings(stats::cor(x[ok], y[ok], method = "spearman"))
  if (is.na(r) || permutations == 0L) return(list(r = r, p = NA_real_))
  mb <- as.matrix(b)
  n <- nrow(mb)
  lower <- lower.tri(mb)
  count <- 0L
  for (k in seq_len(permutations)) {
    perm <- sample.int(n)
    yp <- mb[perm, perm][lower]
    okp <- is.finite(x) & is.finite(yp)
    rp <- suppressWarnings(stats::cor(x[okp], yp[okp], method = "spearman"))
    if (!is.na(rp) && rp >= r) count <- count + 1L
  }
  list(r = r, p = (1 + count) / (1 + permutations))
}

#' Print a comparison of taxonomic sources
#'
#' @param x A `taxodist_source_comparison` returned by [compare_sources()].
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisibly returns `x`.
#' @export
print.taxodist_source_comparison <- function(x, ...) {
  sources <- names(x$resolutions)
  cli::cli_h2("Comparison of taxonomic sources")
  cli::cli_bullets(c(
    "*" = "Sources: {.val {sources}}",
    "*" = "Metric: {.val {x$metric}}",
    "*" = "Taxa resolved in every source: {length(x$shared)} of {nrow(x$taxa)}"
  ))
  depth_cols <- paste0("depth_", sources)
  medians <- vapply(depth_cols, function(col) {
    stats::median(x$taxa[[col]], na.rm = TRUE)
  }, numeric(1))
  cli::cli_text("Median lineage depth: {paste0(sources, ' ', medians, collapse = ', ')}")
  cli::cli_text("")
  out <- x$agreement
  out$spearman <- round(out$spearman, 3)
  print.data.frame(out, row.names = FALSE)
  invisible(x)
}
