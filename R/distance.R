#' Compute the phylogenetic distance between two taxa
#'
#' Given two taxon names, retrieves their lineages from The Taxonomicon and
#' computes a normalized taxonomic distance. The default method is the
#' Jaccard-based distance:
#'
#' \deqn{d_{jaccard}(A, B) = 1 - \frac{\text{depth}(\text{MRCA}(A,B))}
#' {\text{depth}(A) + \text{depth}(B) - \text{depth}(\text{MRCA}(A,B))}}
#'
#' This returns a value in \eqn{[0,1]}, normalized for lineage depth so that
#' deeply resolved taxa such as mammals are not penalized relative to shallowly
#' resolved ones. The underlying raw metric satisfies the triangle inequality
#' \eqn{d(A,C) \leq d(A,B) + d(B,C)}.
#'
#' @param taxon_a A character string giving the first taxon name.
#' @param taxon_b A character string giving the second taxon name.
#' @param method Character. Distance method: `"jaccard"` (default, normalized
#'   0-1), `"norm"` (raw divided by total depth, 0-1), or `"raw"` (integer
#'   tree metric).
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{`distance`}{Numeric. The distance between the two taxa (0-1 for
#'     jaccard/norm, integer for raw).}
#'   \item{`method`}{Character. The method used.}
#'   \item{`mrca`}{Character. The name of the most recent common ancestor.}
#'   \item{`mrca_depth`}{Integer. The depth of the MRCA node.}
#'   \item{`depth_a`}{Integer. The lineage depth of taxon A.}
#'   \item{`depth_b`}{Integer. The lineage depth of taxon B.}
#'   \item{`taxon_a`}{Character. Name of the first taxon.}
#'   \item{`taxon_b`}{Character. Name of the second taxon.}
#' }
#' Returns `NULL` if either taxon cannot be found.
#'
#' @seealso [mrca()], [distance_matrix()], [get_lineage()]
#'
#' @references
#' Brands, S.J. (1989 onwards). Systema Naturae 2000. Amsterdam, The
#' Netherlands. Retrieved from The Taxonomicon,
#' \url{http://taxonomicon.taxonomy.nl}.
#'
#' @export
#' @examples
#' \dontrun{
#' # Distance between two theropods
#' taxo_distance("Tyrannosaurus", "Velociraptor")
#'
#' # Distance between very distantly related taxa
#' taxo_distance("Tyrannosaurus", "Quercus")
#'
#' # Distance between two oviraptorid genera
#' taxo_distance("Nomingia", "Huanansaurus")
#' }
taxo_distance <- function(taxon_a, taxon_b, method = "jaccard", verbose = FALSE) {
  lin_a <- get_lineage(taxon_a, verbose = verbose)
  lin_b <- get_lineage(taxon_b, verbose = verbose)

  if (is.null(lin_a)) {
    cli::cli_alert_danger("Could not retrieve lineage for {taxon_a}")
    return(NULL)
  }
  if (is.null(lin_b)) {
    cli::cli_alert_danger("Could not retrieve lineage for {taxon_b}")
    return(NULL)
  }

  .compute_distance(lin_a, lin_b, taxon_a, taxon_b, method = method)
}

#' Compute the most recent common ancestor of two taxa
#'
#' Retrieves the lineages of two taxa and returns the name of their most
#' recent common ancestor (MRCA) — the deepest node shared by both lineages.
#'
#' @param taxon_a A character string giving the first taxon name.
#' @param taxon_b A character string giving the second taxon name.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A character string giving the name of the MRCA, or `NULL` if
#'   either taxon cannot be found or no common ancestor exists.
#'
#' @export
#' @examples
#' \dontrun{
#' mrca("Tyrannosaurus", "Velociraptor")   # "Coelurosauria"
#' mrca("Tyrannosaurus", "Triceratops")    # "Dinosauria"
#' mrca("Tyrannosaurus", "Homo")           # "Amniota"
#' }
mrca <- function(taxon_a, taxon_b, verbose = FALSE) {
  result <- taxo_distance(taxon_a, taxon_b, verbose = verbose)
  if (is.null(result)) return(NULL)
  result$mrca
}

#' Compute pairwise taxonomic distances for a set of taxa
#'
#' Given a vector of taxon names, computes all pairwise phylogenetic distances
#' and returns a symmetric distance matrix. Lineages are cached after first
#' retrieval to minimise redundant network requests.
#'
#' @param taxa A character vector of taxon names.
#' @param method Character. Distance method passed to [taxo_distance()].
#'   One of `"jaccard"` (default), `"norm"`, or `"raw"`.
#' @param verbose Logical. If `TRUE`, prints progress for each pair.
#'   Default `FALSE`.
#' @param progress Logical. If `TRUE`, shows a progress bar. Default `TRUE`.
#'
#' @return A symmetric numeric matrix of class `"dist"` containing pairwise
#'   distances. Row and column names are set to the input taxon names.
#'   Taxa that could not be found are included with `NA` distances.
#'
#' @seealso [taxo_distance()], [closest_relative()]
#'
#' @export
#' @examples
#' \dontrun{
#' theropods <- c("Tyrannosaurus", "Velociraptor", "Spinosaurus",
#'                "Allosaurus", "Carnotaurus")
#' mat <- distance_matrix(theropods)
#' print(mat)
#' }
distance_matrix <- function(taxa, method = "jaccard", verbose = FALSE, progress = TRUE) {
  n <- length(taxa)
  mat <- matrix(NA_real_, nrow = n, ncol = n,
                dimnames = list(taxa, taxa))
  diag(mat) <- 0

  # pre-fetch all lineages in parallel
  if (progress) cli::cli_alert_info("Fetching {n} lineages in parallel...")
  oplan <- future::plan(future::multisession)
  on.exit(future::plan(oplan), add = TRUE)

  lineages <- furrr::future_map(
    taxa,
    function(t) get_lineage(t, verbose = verbose),
    .options = furrr::furrr_options(seed = NULL)
  )
  names(lineages) <- taxa
  if (progress) cli::cli_alert_success("Lineages fetched.")

  # compute pairwise distances
  total_pairs <- n * (n - 1) / 2
  if (progress) cli::cli_progress_bar("Computing distances", total = total_pairs)

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      if (!is.null(lineages[[i]]) && !is.null(lineages[[j]])) {
        result <- .compute_distance(
          lineages[[i]], lineages[[j]], taxa[i], taxa[j], method = method
        )
        mat[i, j] <- mat[j, i] <- result$distance
      }
      if (progress) cli::cli_progress_update()
    }
  }
  if (progress) cli::cli_progress_done()

  stats::as.dist(mat)
}

#' Find the closest relative of a taxon among a set of candidates
#'
#' Given a query taxon and a vector of candidate taxa, returns the candidate
#' with the smallest phylogenetic distance to the query.
#'
#' @param taxon A character string giving the query taxon name.
#' @param candidates A character vector of candidate taxon names to compare
#'   against.
#' @param method Character. Distance method passed to [taxo_distance()].
#'   One of `"jaccard"` (default), `"norm"`, or `"raw"`.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A data frame with columns `taxon` (candidate name) and `distance`
#'   (tree metric distance), sorted by distance ascending. Returns `NULL` if
#'   the query taxon cannot be found.
#'
#' @export
#' @examples
#' \dontrun{
#' closest_relative("Tyrannosaurus",
#'   c("Velociraptor", "Triceratops", "Brachiosaurus", "Allosaurus"))
#' }
closest_relative <- function(taxon, candidates, method = "jaccard", verbose = FALSE) {
  query_lin <- get_lineage(taxon, verbose = verbose)
  if (is.null(query_lin)) {
    cli::cli_alert_danger("Could not retrieve lineage for {taxon}")
    return(NULL)
  }

  oplan <- future::plan(future::multisession)
  on.exit(future::plan(oplan), add = TRUE)

  results <- furrr::future_map_dfr(
    candidates,
    function(cand) {
      cand_lin <- get_lineage(cand, verbose = verbose)
      if (is.null(cand_lin)) {
        return(data.frame(taxon = cand, distance = NA_real_))
      }
      dist_result <- .compute_distance(query_lin, cand_lin, taxon, cand,
                                       method = method)
      data.frame(taxon = cand, distance = dist_result$distance)
    },
    .options = furrr::furrr_options(seed = NULL)
  )

  results[order(results$distance, na.last = TRUE), ]
}

#' Get the lineage depth of a taxon
#'
#' Returns the number of nodes in the lineage of a taxon, from root to tip.
#' This reflects how deeply nested the taxon is within the taxonomic hierarchy.
#'
#' @param taxon A character string giving the taxon name.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return An integer giving the lineage depth, or `NULL` if the taxon cannot
#'   be found.
#'
#' @export
#' @examples
#' \dontrun{
#' lineage_depth("Tyrannosaurus")   # deep — many intermediate clades
#' lineage_depth("Biota")           # shallow — near root
#' }
lineage_depth <- function(taxon, verbose = FALSE) {
  lin <- get_lineage(taxon, verbose = verbose)
  if (is.null(lin)) return(NULL)
  length(lin)
}

#' Check whether a taxon is covered by The Taxonomicon
#'
#' Queries The Taxonomicon for a taxon name and returns a logical indicating
#' whether the taxon was found. Useful for pre-screening a list of names
#' before running distance computations.
#'
#' @param taxa A character vector of one or more taxon names.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A named logical vector. `TRUE` indicates the taxon was found,
#'   `FALSE` indicates it was not.
#'
#' @export
#' @examples
#' \dontrun{
#' check_coverage(c("Tyrannosaurus", "Velociraptor", "Fakeosaurus"))
#' }
check_coverage <- function(taxa, verbose = FALSE) {
  result <- purrr::map_lgl(taxa, function(t) {
    !is.null(get_taxonomicon_id(t, verbose = verbose))
  })
  names(result) <- taxa
  result
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.compute_distance <- function(lin_a, lin_b, name_a = "A", name_b = "B",
                               method = "jaccard") {
  # find longest common prefix (MRCA)
  min_len    <- min(length(lin_a), length(lin_b))
  mrca_depth <- 0L

  for (i in seq_len(min_len)) {
    if (lin_a[i] == lin_b[i]) {
      mrca_depth <- i
    } else {
      break
    }
  }

  depth_a   <- length(lin_a)
  depth_b   <- length(lin_b)
  mrca_name <- if (mrca_depth > 0L) lin_a[mrca_depth] else NA_character_
  d_raw     <- depth_a + depth_b - 2L * mrca_depth

  distance <- switch(method,
    raw     = d_raw,
    norm    = d_raw / (depth_a + depth_b),
    jaccard = 1 - mrca_depth / (depth_a + depth_b - mrca_depth),
    stop("method must be 'raw', 'norm', or 'jaccard'")
  )

  list(
    distance   = distance,
    method     = method,
    mrca       = mrca_name,
    mrca_depth = mrca_depth,
    depth_a    = depth_a,
    depth_b    = depth_b,
    taxon_a    = name_a,
    taxon_b    = name_b
  )
}
