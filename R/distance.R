#' Compute the phylogenetic distance between two taxa
#'
#' Given two taxon names, retrieves their lineages from The Taxonomicon and
#' computes a taxonomic distance based on the depth of their most recent
#' common ancestor (MRCA):
#'
#' \deqn{d(A, B) = \frac{1}{\text{depth}(\text{MRCA}(A,B))}}
#'
#' The deeper the shared ancestor, the smaller (closer to zero) the distance.
#' This metric ensures that taxa diverging at the same node are always
#' equidistant from any third taxon, regardless of lineage depth differences
#' below the split.
#'
#' @param taxon_a A character string giving the first taxon name.
#' @param taxon_b A character string giving the second taxon name.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A named list of class `"taxodist_result"` with the following elements:
#' \describe{
#'   \item{`distance`}{Numeric. The distance between the two taxa. Returns 0
#'     if one taxon is an ancestor of the other.}
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
#' \donttest{
#' # Distance between two theropods
#' taxo_distance("Tyrannosaurus", "Velociraptor")
#'
#' # Distance between very distantly related taxa
#' taxo_distance("Tyrannosaurus", "Quercus")
#'
#' # Distance between two oviraptorid genera
#' taxo_distance("Nomingia", "Huanansaurus")
#' }
taxo_distance <- function(taxon_a, taxon_b, verbose = FALSE) {
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

  .compute_distance(lin_a, lin_b, taxon_a, taxon_b)
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
#' \donttest{
#' mrca("Tyrannosaurus", "Velociraptor")   # "Tyrannoraptora"
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
#' \donttest{
#' theropods <- c("Tyrannosaurus", "Velociraptor", "Spinosaurus",
#'                "Allosaurus", "Carnotaurus")
#' mat <- distance_matrix(theropods)
#' print(mat)
#' }
distance_matrix <- function(taxa, verbose = FALSE, progress = TRUE) {
  n <- length(taxa)
  mat <- matrix(NA_real_, nrow = n, ncol = n,
                dimnames = list(taxa, taxa))
  diag(mat) <- 0

  # fetch lineages sequentially in main process (cache is shared)
  if (progress) cli::cli_alert_info("Fetching {n} lineages...")
  lineages <- lapply(taxa, function(t) get_lineage(t, verbose = verbose))
  names(lineages) <- taxa
  if (progress) cli::cli_alert_success("Lineages fetched.")

  # compute pairwise distances
  total_pairs <- n * (n - 1) / 2
  if (progress) cli::cli_progress_bar("Computing distances", total = total_pairs)
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      if (!is.null(lineages[[i]]) && !is.null(lineages[[j]])) {
        result <- .compute_distance(
          lineages[[i]], lineages[[j]], taxa[i], taxa[j]
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
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A data frame with columns `taxon` (candidate name) and `distance`
#'   (tree metric distance), sorted by distance ascending. Returns `NULL` if
#'   the query taxon cannot be found.
#'
#' @export
#' @examples
#' \donttest{
#' closest_relative("Tyrannosaurus",
#'   c("Velociraptor", "Triceratops", "Brachiosaurus", "Allosaurus"))
#' }
closest_relative <- function(taxon, candidates, verbose = FALSE) {
  query_lin <- get_lineage(taxon, verbose = verbose)
  if (is.null(query_lin)) {
    cli::cli_alert_danger("Could not retrieve lineage for {taxon}")
    return(NULL)
  }

  results <- do.call(rbind, purrr::map(
    candidates,
    function(cand) {
      cand_lin <- get_lineage(cand, verbose = verbose)
      if (is.null(cand_lin)) {
        return(data.frame(taxon = cand, distance = NA_real_))
      }
      dist_result <- .compute_distance(query_lin, cand_lin, taxon, cand)
      data.frame(taxon = cand, distance = dist_result$distance)
    }
  ))

  results[order(results$distance, na.last = TRUE), ]
}

#' Compute distances from a focal taxon to a community of taxa
#'
#' Given a focal taxon and a vector of community taxa, retrieves all lineages
#' and computes the taxonomic distance from the focal taxon to each member of
#' the community. Returns a sorted data frame that also reports the most recent
#' common ancestor (MRCA) for each pair, making it easy to interpret *why* a
#' taxon is close or distant.
#'
#' @param focal A character string giving the focal taxon name.
#' @param community A character vector of community taxon names to compare
#'   against. The focal taxon may be included — it will receive a distance
#'   of `0`.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#' @param progress Logical. If `TRUE`, shows a progress bar. Default `TRUE`.
#'
#' @return A data frame of class `"taxodist_focal"` with columns:
#' \describe{
#'   \item{`taxon`}{Character. Community taxon name.}
#'   \item{`distance`}{Numeric. Taxonomic distance to the focal taxon.}
#'   \item{`mrca`}{Character. Name of the most recent common ancestor.}
#'   \item{`mrca_depth`}{Integer. Depth of the MRCA node.}
#' }
#' Rows are sorted by `distance` ascending (closest relatives first).
#' Returns `NULL` if the focal taxon cannot be found.
#'
#' @seealso [closest_relative()], [distance_matrix()], [taxo_distance()]
#' @export
#' @examples
#' \donttest{
#' community <- c("Velociraptor", "Triceratops", "Brachiosaurus",
#'                "Allosaurus", "Homo", "Quercus")
#' focal_distances("Tyrannosaurus", community)
#' }
focal_distances <- function(focal, community, verbose = FALSE, progress = TRUE) {
  focal_lin <- get_lineage(focal, verbose = verbose)
  if (is.null(focal_lin)) {
    cli::cli_alert_danger("Could not retrieve lineage for {focal}")
    return(NULL)
  }

  if (progress) cli::cli_progress_bar("Computing focal distances", total = length(community))

  results <- vector("list", length(community))
  for (i in seq_along(community)) {
    taxon <- community[[i]]

    if (identical(taxon, focal)) {
      results[[i]] <- data.frame(
        taxon      = taxon,
        distance   = 0,
        mrca       = focal,
        mrca_depth = length(focal_lin),
        stringsAsFactors = FALSE
      )
    } else {
      cand_lin <- get_lineage(taxon, verbose = verbose)
      if (is.null(cand_lin)) {
        results[[i]] <- data.frame(
          taxon      = taxon,
          distance   = NA_real_,
          mrca       = NA_character_,
          mrca_depth = NA_integer_,
          stringsAsFactors = FALSE
        )
      } else {
        res <- .compute_distance(focal_lin, cand_lin, focal, taxon)
        results[[i]] <- data.frame(
          taxon      = taxon,
          distance   = res$distance,
          mrca       = res$mrca,
          mrca_depth = res$mrca_depth,
          stringsAsFactors = FALSE
        )
      }
    }

    if (progress) cli::cli_progress_update()
  }

  if (progress) cli::cli_progress_done()

  out <- do.call(rbind, results)
  out <- out[order(out$distance, na.last = TRUE), ]
  rownames(out) <- NULL

  structure(out, class = c("taxodist_focal", "data.frame"), focal = focal)
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
#' \donttest{
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
#' \donttest{
#' check_coverage(c("Tyrannosaurus", "Velociraptor", "Fakeosaurus"))
#' }
check_coverage <- function(taxa, verbose = FALSE) {
  result <- purrr::map_lgl(taxa, function(t) {
    !is.null(get_taxonomicon_id(t, verbose = verbose))
  })
  names(result) <- taxa
  result
}

#' Cluster taxa by taxonomic distance
#'
#' Computes pairwise taxonomic distances and performs hierarchical clustering.
#'
#' @param taxa A character vector of taxon names, or a `dist` object from
#'[distance_matrix()].
#' @param method Clustering method passed to [stats::hclust()]. Default
#'   `"average"` (UPGMA), which works well with taxonomic distances.
#' @param ... Additional arguments passed to [distance_matrix()] (e.g.
#'   `verbose`, `progress`).
#'
#' @return An object of class `"taxodist_cluster"` — a list with:
#' \describe{
#'   \item{`hclust`}{The [stats::hclust()] result.}
#'   \item{`dist`}{The underlying distance matrix.}
#' }
#'
#' @seealso[taxo_ordinate()], [distance_matrix()]
#' @export
#' @examples
#' \donttest{
#' taxa <- c("Tyrannosaurus", "Velociraptor", "Triceratops", "Brachiosaurus")
#' cl <- taxo_cluster(taxa)
#' if (!is.null(cl$hclust)) {
#'   plot(cl)
#' }
#' }
taxo_cluster <- function(taxa, method = "average", ...) {
  d <- if (inherits(taxa, "dist")) taxa else distance_matrix(taxa, ...)
  if (any(is.na(d))) {
    cli::cli_warn("Distance matrix contains NA values (taxa not found or server offline). Clustering skipped.")
    return(structure(list(hclust = NULL, dist = d), class = "taxodist_cluster"))
  }
  hc <- stats::hclust(d, method = method)
  structure(list(hclust = hc, dist = d), class = "taxodist_cluster")
}

#' Ordinate taxa in taxonomic distance space
#'
#' Computes pairwise taxonomic distances and applies classical multidimensional
#' scaling (PCoA) to project taxa into a low-dimensional space.
#'
#' @param taxa A character vector of taxon names, or a `dist` object from
#'[distance_matrix()].
#' @param k Number of dimensions. Default `2`.
#' @param ... Additional arguments passed to [distance_matrix()].
#'
#' @return An object of class `"taxodist_ord"` — a list with:
#' \describe{
#'   \item{`points`}{A matrix of coordinates (taxa x k dimensions).}
#'   \item{`dist`}{The underlying distance matrix.}
#'   \item{`GOF`}{Goodness-of-fit from[stats::cmdscale()].}
#'   \item{`eig`}{The eigenvalues computed during PCoA.}
#' }
#'
#' @seealso [taxo_cluster()], [distance_matrix()]
#' @export
#' @examples
#' \donttest{
#' taxa <- c("Tyrannosaurus", "Velociraptor", "Triceratops", "Brachiosaurus")
#' ord <- taxo_ordinate(taxa)
#' if (!is.null(ord$points)) {
#'   plot(ord$points, type = "n")
#'   text(ord$points, labels = rownames(ord$points))
#' }
#' }
taxo_ordinate <- function(taxa, k = 2, ...) {
  d <- if (inherits(taxa, "dist")) taxa else distance_matrix(taxa, ...)
  if (any(is.na(d))) {
    cli::cli_warn("Distance matrix contains NA values. Ordination skipped.")
    return(structure(list(points = NULL, dist = d, GOF = NULL, eig = NULL), class = "taxodist_ord"))
  }
  cmd <- stats::cmdscale(d, k = k, eig = TRUE)
  structure(list(
    points = cmd$points,
    dist   = d,
    GOF    = cmd$GOF,
    eig    = cmd$eig
  ), class = "taxodist_ord")
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.compute_distance <- function(lin_a, lin_b, name_a = "A", name_b = "B") {
  depth_a <- length(lin_a)
  depth_b <- length(lin_b)

  # find deepest shared node (MRCA) using set intersection
  shared <- intersect(lin_a, lin_b)

  if (length(shared) == 0) {
    return(structure(list(
      distance   = Inf,
      mrca       = NA_character_,
      mrca_depth = 0L,
      depth_a    = depth_a,
      depth_b    = depth_b,
      taxon_a    = name_a,
      taxon_b    = name_b
    ), class = "taxodist_result"))
  }

  # find position of each shared node in lin_a, take the deepest
  positions_in_a <- match(shared, lin_a)
  mrca_idx       <- which.max(positions_in_a)
  mrca_depth     <- positions_in_a[mrca_idx]
  mrca_name      <- lin_a[mrca_depth]
  is_ancestral   <- (mrca_name == lin_a[depth_a]) || (mrca_name == lin_b[depth_b])
  distance       <- if (is_ancestral) 0 else 1 / mrca_depth

  structure(list(
    distance   = distance,
    mrca       = mrca_name,
    mrca_depth = mrca_depth,
    depth_a    = depth_a,
    depth_b    = depth_b,
    taxon_a    = name_a,
    taxon_b    = name_b
  ), class = "taxodist_result")
}
