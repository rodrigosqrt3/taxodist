#' Compute the taxonomic hierarchy distance between two taxa
#'
#' Given two taxon names, retrieves their lineages from The Taxonomicon and
#' computes a taxonomic distance based on the depth of their most recent
#' common ancestor (MRCA):
#'
#' \deqn{d(A,A) = 0}{d(A,A) = 0}
#' \deqn{d(A,B) = 1 / h(A,B)}{d(A,B) = 1 / h(A,B)}
#'
#' Here, `h(A,B)` is the depth of the MRCA of `A` and `B`. The deeper the
#' shared ancestor, the smaller (closer to zero) the distance.
#' Distinct nodes always have positive distance, including when one taxon is
#' an ancestor of the other. This makes the measure an ultrametric on each
#' connected taxonomic hierarchy. Use [is_member()] when the question is
#' whether one taxon belongs to a clade rather than how far apart the two
#' hierarchy nodes are.
#'
#' @param taxon_a A character string giving the first taxon name.
#' @param taxon_b A character string giving the second taxon name.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A named list of class `"taxodist_result"` with the following elements:
#' \describe{
#'   \item{`distance`}{Numeric. The distance between the two taxa. Returns 0
#'     only when their complete lineages are identical.}
#'   \item{`mrca`}{Character. The name of the most recent common ancestor.}
#'   \item{`mrca_depth`}{Integer. The depth of the MRCA node.}
#'   \item{`depth_a`}{Integer. The lineage depth of taxon A.}
#'   \item{`depth_b`}{Integer. The lineage depth of taxon B.}
#'   \item{`taxon_a`}{Character. Name of the first taxon.}
#'   \item{`taxon_b`}{Character. Name of the second taxon.}
#' }
#' Returns `NULL` if either taxon cannot be found.
#'
#' @seealso [mrca()], [distance_matrix()], [get_lineage()], [is_member()],
#'   [taxo_path()]
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
#' Given a vector of taxon names, computes all pairwise taxonomic hierarchy distances
#' and returns a symmetric distance matrix. Lineages are cached after first
#' retrieval to minimise redundant network requests.
#'
#' @param taxa A character vector of taxon names, a `taxodist_resolution`
#'   object returned by [taxo_resolve()], or a `taxodist_bundle`. Resolution
#'   objects reuse their stored lineages and retain the original inputs as
#'   matrix labels; bundles return their stored matrix.
#' @param verbose Logical. If `TRUE`, prints progress for each pair.
#'   Default `FALSE`.
#' @param progress Logical. If `TRUE`, shows a progress bar. Default `TRUE`.
#'
#' @return A symmetric numeric matrix of class `"dist"` containing pairwise
#'   distances. Row and column names are set to the input taxon names.
#'   Taxa that could not be found are included with `NA` distances.
#'   An empty input returns an empty `"dist"` object.
#'
#' @details
#' Distances are computed in vectorized form: every lineage node is first
#' replaced by an integer identifying its complete path from the root, so the
#' MRCA depth of two taxa is simply the number of depths at which these
#' identifiers agree. This gives exactly the same values as
#' [taxo_distance()] while scaling to thousands of taxa.
#'
#' @seealso [taxo_distance()], [cross_distances()], [closest_relative()],
#'   [taxo_resolve()], [taxo_bundle()], [taxo_tree()]
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
  if (inherits(taxa, "taxodist_bundle")) {
    .validate_taxodist_bundle(taxa, check_distances = FALSE)
    return(taxa$matrix)
  }

  input <- .matrix_input(taxa, get_lineage, verbose = verbose, progress = progress)
  labels <- input$labels
  n <- length(labels)

  values <- if (n > 1L) {
    .pairwise_distances(input$lineages, progress = progress)
  } else {
    numeric(0)
  }

  structure(
    values,
    Size = n,
    Labels = labels,
    Diag = FALSE,
    Upper = FALSE,
    class = "dist"
  )
}

#' Compute taxonomic distances between two sets of taxa
#'
#' Computes the rectangular matrix of distances from every taxon in `x` to
#' every taxon in `y`, for example from the species of a community to a set
#' of reference taxa. This avoids computing the full square matrix of the
#' combined set.
#'
#' @param x,y Character vectors of taxon names, or `taxodist_resolution`
#'   objects returned by [taxo_resolve()] or [taxo_from_lineages()].
#' @param verbose Logical. If `TRUE`, prints retrieval messages.
#' @param progress Logical. If `TRUE`, shows progress messages.
#'
#' @return A numeric matrix with one row per element of `x` and one column
#'   per element of `y`, labeled with the inputs. Taxa that could not be
#'   found have `NA` distances; taxa without a common ancestor have `Inf`.
#'
#' @seealso [distance_matrix()], [focal_distances()]
#' @export
#' @examples
#' community <- taxo_from_lineages(list(
#'   Alpha = c("Biota", "Animalia", "Chordata", "Alpha"),
#'   Beta  = c("Biota", "Animalia", "Arthropoda", "Beta")
#' ))
#' reference <- taxo_from_lineages(list(
#'   Gamma = c("Biota", "Animalia", "Chordata", "Gamma"),
#'   Delta = c("Biota", "Plantae", "Delta")
#' ))
#' cross_distances(community, reference)
cross_distances <- function(x, y, verbose = FALSE, progress = TRUE) {
  input_x <- .matrix_input(
    x, get_lineage, verbose = verbose, progress = progress, arg = "x"
  )
  input_y <- .matrix_input(
    y, get_lineage, verbose = verbose, progress = progress, arg = "y"
  )
  nx <- length(input_x$labels)
  ny <- length(input_y$labels)

  out <- matrix(
    NA_real_, nrow = nx, ncol = ny,
    dimnames = list(input_x$labels, input_y$labels)
  )
  if (nx == 0L || ny == 0L) return(out)

  lineages <- c(input_x$lineages, input_y$lineages)
  ids <- .prefix_node_ids(lineages)
  lens <- lengths(lineages)
  missing <- vapply(lineages, is.null, logical(1))
  x_idx <- seq_len(nx)

  for (j in seq_len(ny)) {
    col <- nx + j
    if (missing[[col]]) next
    out[, j] <- .distances_to(ids, lens, x_idx, col)
  }
  out[missing[x_idx], ] <- NA_real_
  out
}

#' Find the closest relative of a taxon among a set of candidates
#'
#' Given a query taxon and a vector of candidate taxa, returns the candidate
#' with the smallest taxonomic hierarchy distance to the query.
#'
#' @param taxon A character string giving the query taxon name.
#' @param candidates A character vector of candidate taxon names to compare
#'   against.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A data frame with columns `taxon` (candidate name) and `distance`
#'   (ultrametric distance), sorted by distance ascending. Returns `NULL` if
#'   the query taxon cannot be found. An empty candidate vector returns an
#'   empty data frame with the same columns.
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

  distances <- vapply(candidates, function(cand) {
    cand_lin <- get_lineage(cand, verbose = verbose)
    if (is.null(cand_lin)) return(NA_real_)
    .compute_distance(query_lin, cand_lin, taxon, cand)$distance
  }, numeric(1), USE.NAMES = FALSE)

  results <- data.frame(
    taxon = as.character(candidates),
    distance = distances,
    stringsAsFactors = FALSE
  )
  results <- results[order(results$distance, na.last = TRUE), , drop = FALSE]
  rownames(results) <- NULL
  results
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
#' Returns `NULL` if the focal taxon cannot be found. An empty community
#' returns an empty object with the same columns and class.
#'
#' @seealso [closest_relative()], [cross_distances()], [distance_matrix()],
#'   [taxo_distance()]
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

  if (length(community) == 0L) {
    out <- data.frame(
      taxon = character(0),
      distance = numeric(0),
      mrca = character(0),
      mrca_depth = integer(0),
      stringsAsFactors = FALSE
    )
    return(structure(
      out,
      class = c("taxodist_focal", "data.frame"),
      focal = focal
    ))
  }

  if (progress) cli::cli_progress_bar("Computing focal distances", total = length(community))

  results <- vector("list", length(community))
  for (i in seq_along(community)) {
    taxon <- community[[i]]
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
#'   `FALSE` that it is absent from The Taxonomicon, and `NA` that the
#'   answer is unknown because the server could not be reached.
#'
#' @seealso [taxo_resolve()] for a complete audit of each name.
#' @export
#' @examples
#' \donttest{
#' check_coverage(c("Tyrannosaurus", "Velociraptor", "Fakeosaurus"))
#' }
check_coverage <- function(taxa, verbose = FALSE) {
  result <- vapply(taxa, function(t) {
    lookup <- .lookup_taxonomicon_id(t, verbose = verbose)
    switch(lookup$status,
      found = TRUE,
      not_found = FALSE,
      NA
    )
  }, logical(1), USE.NAMES = FALSE)
  names(result) <- taxa
  result
}

#' Cluster taxa by taxonomic distance
#'
#' Computes pairwise taxonomic distances and performs hierarchical clustering.
#'
#' @param taxa A character vector of taxon names, or a `dist` object from
#'   [distance_matrix()].
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
#' If fewer than two taxa are supplied, `hclust` is `NULL` and the distance
#' object is preserved.
#'
#' @seealso [taxo_tree()] for a tree whose internal nodes are the named clades
#'   of the classification, [taxo_ordinate()], [distance_matrix()]
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
  if (any(!is.finite(d))) {
    cli::cli_warn("Distance matrix contains infinite values (no shared ancestor). Clustering skipped.")
    return(structure(list(hclust = NULL, dist = d), class = "taxodist_cluster"))
  }
  if (attr(d, "Size") < 2L) {
    cli::cli_warn("At least two taxa are required for clustering. Clustering skipped.")
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
#'   [distance_matrix()].
#' @param k Number of dimensions. Default `2`.
#' @param ... Additional arguments passed to [distance_matrix()].
#'
#' @return An object of class `"taxodist_ord"` — a list with:
#' \describe{
#'   \item{`points`}{A matrix of coordinates (taxa x k dimensions).}
#'   \item{`dist`}{The underlying distance matrix.}
#'   \item{`GOF`}{Goodness-of-fit from [stats::cmdscale()].}
#'   \item{`eig`}{The eigenvalues computed during PCoA.}
#' }
#' At least two taxa are required. When `k` is greater than the maximum
#' possible dimension (`n - 1`), it is reduced automatically.
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
  if (any(!is.finite(d))) {
    cli::cli_warn("Distance matrix contains infinite values (no shared ancestor). Ordination skipped.")
    return(structure(list(points = NULL, dist = d, GOF = NULL, eig = NULL), class = "taxodist_ord"))
  }
  n_taxa <- attr(d, "Size")
  if (n_taxa < 2L) {
    cli::cli_warn("At least two taxa are required for ordination. Ordination skipped.")
    return(structure(list(points = NULL, dist = d, GOF = NULL, eig = NULL), class = "taxodist_ord"))
  }
  valid_k <- length(k) == 1L &&
    (is.integer(k) || is.double(k)) &&
    is.finite(k) && k >= 1 && k == floor(k)
  if (!valid_k) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }
  max_k <- n_taxa - 1L
  if (k > max_k) {
    cli::cli_warn("`k` reduced from {k} to {max_k}, the maximum for {n_taxa} taxa.")
    k <- max_k
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

# Turn the supported inputs (character vector, taxodist_resolution, or
# taxodist_bundle) into parallel vectors of labels and lineages. Lineages of
# character inputs are retrieved with `fetch` (get_lineage()) and cached under
# their normalized name. Callers pass `fetch` explicitly so that it is
# resolved in their own environment.
#' @keywords internal
.matrix_input <- function(taxa, fetch, verbose = FALSE, progress = TRUE,
                          arg = "taxa") {
  if (inherits(taxa, "taxodist_bundle")) {
    taxa <- taxa$resolution
  }
  if (inherits(taxa, "taxodist_resolution")) {
    required <- c("input", "status", "lineage")
    if (!all(required %in% names(taxa))) {
      cli::cli_abort("Invalid {.cls taxodist_resolution} object.")
    }
    return(list(labels = taxa$input, lineages = unclass(taxa$lineage)))
  }
  if (is.list(taxa)) {
    cli::cli_abort(
      "{.arg {arg}} must be a character vector or a {.cls taxodist_resolution}."
    )
  }
  taxa <- as.character(taxa)
  if (length(taxa) == 0L) {
    return(list(labels = taxa, lineages = list()))
  }

  # Fetch lineages sequentially in the main process (cache is shared).
  if (progress) cli::cli_alert_info("Fetching {length(taxa)} lineages...")
  lineages <- lapply(taxa, function(taxon) {
    if (is.na(taxon)) return(NULL)
    cache_key <- paste0("matrix_lineage_", .normalize_taxon(taxon))
    if (exists(cache_key, envir = .taxodist_cache, inherits = FALSE)) {
      return(get(cache_key, envir = .taxodist_cache, inherits = FALSE))
    }

    lineage <- fetch(taxon, verbose = verbose)
    if (!is.null(lineage)) {
      assign(cache_key, lineage, envir = .taxodist_cache)
    }
    lineage
  })
  if (progress) cli::cli_alert_success("Lineages fetched.")
  list(labels = taxa, lineages = lineages)
}

# Replace each lineage node by an integer that identifies its complete path
# from the root. Two lineages share the node at depth d if and only if their
# identifiers at depth d are equal, and equality at depth d implies equality
# at every shallower depth, so the MRCA depth of two taxa equals the number
# of depths at which their identifiers agree. Missing depths are NA.
#' @keywords internal
.prefix_node_ids <- function(lineages) {
  n <- length(lineages)
  lens <- lengths(lineages)
  max_depth <- if (n > 0L) max(lens) else 0L
  ids <- matrix(NA_integer_, nrow = n, ncol = max_depth)
  if (max_depth == 0L) return(ids)

  flat <- unlist(lineages, use.names = FALSE)
  owner <- rep.int(seq_len(n), lens)
  depth <- sequence(lens)
  parent <- integer(n)
  for (d in seq_len(max_depth)) {
    at_depth <- depth == d
    rows <- owner[at_depth]
    key <- paste(parent[rows], flat[at_depth], sep = "\r")
    node <- match(key, unique(key))
    ids[rows, d] <- node
    parent[rows] <- node
  }
  ids
}

# Distances from the lineage in row `j` of `ids` to the lineages in rows `i`.
#' @keywords internal
.distances_to <- function(ids, lens, i, j) {
  shared <- ids[i, , drop = FALSE] == rep(ids[j, ], each = length(i))
  mrca_depth <- rowSums(shared, na.rm = TRUE)
  values <- 1 / mrca_depth
  identical_lineage <- mrca_depth > 0 &
    mrca_depth == lens[[j]] & lens[i] == lens[[j]]
  values[identical_lineage] <- 0
  values
}

# Lower-triangle distances in base `dist` storage order (column by column).
# NULL lineages (taxa that could not be resolved) yield NA.
#' @keywords internal
.pairwise_distances <- function(lineages, progress = FALSE) {
  n <- length(lineages)
  ids <- .prefix_node_ids(lineages)
  lens <- lengths(lineages)
  missing <- vapply(lineages, is.null, logical(1))

  values <- numeric(n * (n - 1) / 2)
  position <- 0L
  if (progress) {
    cli::cli_progress_bar("Computing distances", total = n - 1L)
  }
  for (column in seq_len(n - 1L)) {
    rows <- (column + 1L):n
    block <- .distances_to(ids, lens, rows, column)
    block[missing[rows] | missing[[column]]] <- NA_real_
    values[position + seq_along(rows)] <- block
    position <- position + length(rows)
    if (progress) cli::cli_progress_update()
  }
  if (progress) cli::cli_progress_done()
  values
}

#' @keywords internal
.common_prefix_depth <- function(lin_a, lin_b) {
  common_length <- min(length(lin_a), length(lin_b))
  if (common_length == 0L) return(0L)

  common_idx <- seq_len(common_length)
  matches <- lin_a[common_idx] == lin_b[common_idx]
  first_mismatch <- which(is.na(matches) | !matches)

  if (length(first_mismatch) > 0L) {
    first_mismatch[[1L]] - 1L
  } else {
    common_length
  }
}

#' @keywords internal
.compute_distance <- function(lin_a, lin_b, name_a = "A", name_b = "B") {
  depth_a <- length(lin_a)
  depth_b <- length(lin_b)

  # The shared ancestry must be a continuous prefix from the root. A taxon
  # name repeated after the lineages diverge is a homonym, not a shared node.
  mrca_depth <- .common_prefix_depth(lin_a, lin_b)

  if (mrca_depth == 0L) {
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

  mrca_name    <- lin_a[mrca_depth]
  same_lineage <- depth_a == depth_b && mrca_depth == depth_a
  distance     <- if (same_lineage) 0 else 1 / mrca_depth

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
