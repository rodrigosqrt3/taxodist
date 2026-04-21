#' Print method for taxodist distance results
#'
#' @param x A list returned by [taxo_distance()].
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{x}. Called for side effects (printing).
#' @export
print.taxodist_result <- function(x, ...) {
  cli::cli_h2("Taxonomic Distance")
  cli::cli_bullets(c(
    "*" = "{x$taxon_a} vs {x$taxon_b}",
    " " = "Distance : {x$distance}",
    " " = "MRCA     : {x$mrca} (depth {x$mrca_depth})",
    " " = "Depth A  : {x$depth_a}",
    " " = "Depth B  : {x$depth_b}"
  ))
  invisible(x)
}

#' Compare lineages of two taxa side by side
#'
#' Prints the lineages of two taxa aligned at their most recent common
#' ancestor, making the point of divergence easy to identify.
#'
#' @param taxon_a A character string giving the first taxon name.
#' @param taxon_b A character string giving the second taxon name.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return Invisibly returns a list with elements `lineage_a`, `lineage_b`,
#'   and `mrca_depth`.
#'
#' @export
#' @examples
#' \donttest{
#' compare_lineages("Tyrannosaurus", "Velociraptor")
#' compare_lineages("Tyrannosaurus", "Triceratops")
#' }
compare_lineages <- function(taxon_a, taxon_b, verbose = FALSE) {
  lin_a <- get_lineage(taxon_a, verbose = verbose)
  lin_b <- get_lineage(taxon_b, verbose = verbose)

  if (is.null(lin_a) || is.null(lin_b)) {
    cli::cli_alert_danger("Could not retrieve one or both lineages")
    return(invisible(NULL))
  }

  result <- .compute_distance(lin_a, lin_b, taxon_a, taxon_b)
  mrca_d <- result$mrca_depth

  cli::cli_h2("Lineage Comparison")
  cli::cli_text("MRCA: {result$mrca} at depth {mrca_d}")
  cli::cli_text("")

  # shared trunk
  if (mrca_d > 0) {
    cli::cli_text("{.strong Shared lineage ({mrca_d} nodes):}")
    for (node in lin_a[seq_len(mrca_d)]) {
      cli::cli_bullets(c(" " = "{node}"))
    }
  }

  # divergence
  cli::cli_text("")
  cli::cli_text("{.strong {taxon_a} only ({length(lin_a) - mrca_d} nodes):}")
  if (length(lin_a) > mrca_d) {
    for (node in lin_a[(mrca_d + 1):length(lin_a)]) {
      cli::cli_bullets(c(">" = "{node}"))
    }
  }

  cli::cli_text("")
  cli::cli_text("{.strong {taxon_b} only ({length(lin_b) - mrca_d} nodes):}")
  if (length(lin_b) > mrca_d) {
    for (node in lin_b[(mrca_d + 1):length(lin_b)]) {
      cli::cli_bullets(c(">" = "{node}"))
    }
  }

  invisible(list(
    lineage_a  = lin_a,
    lineage_b  = lin_b,
    mrca_depth = mrca_d
  ))
}

#' List all clades shared between two taxa
#'
#' Returns the vector of clade names forming the shared trunk of two taxa's
#' lineages, from root down to (and including) their MRCA.
#'
#' @param taxon_a A character string giving the first taxon name.
#' @param taxon_b A character string giving the second taxon name.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A character vector of shared clade names ordered from root to MRCA,
#'   or `NULL` if either taxon cannot be found.
#'
#' @export
#' @examples
#' \donttest{
#' shared_clades("Tyrannosaurus", "Velociraptor")
#' shared_clades("Tyrannosaurus", "Homo")
#' }
shared_clades <- function(taxon_a, taxon_b, verbose = FALSE) {
  lin_a <- get_lineage(taxon_a, verbose = verbose)
  lin_b <- get_lineage(taxon_b, verbose = verbose)
  if (is.null(lin_a) || is.null(lin_b)) return(NULL)

  result <- .compute_distance(lin_a, lin_b, taxon_a, taxon_b)
  if (result$mrca_depth == 0) return(character(0))
  lin_a[seq_len(result$mrca_depth)]
}

#' Test whether one taxon is nested within another
#'
#' Returns `TRUE` if `taxon` is a member of `clade` — i.e., if the clade name
#' appears in the taxon's lineage.
#'
#' @param taxon A character string giving the taxon name to test.
#' @param clade A character string giving the clade name to test membership in.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A logical value, or `NULL` if the taxon cannot be found.
#'
#' @export
#' @examples
#' \donttest{
#' is_member("Tyrannosaurus", "Theropoda")    # TRUE
#' is_member("Triceratops", "Theropoda")      # FALSE
#' is_member("Homo", "Amniota")               # TRUE
#' }
is_member <- function(taxon, clade, verbose = FALSE) {
  lin <- get_lineage(taxon, verbose = verbose)
  if (is.null(lin)) return(NULL)
  any(grepl(paste0("^", clade), lin, ignore.case = TRUE))
}

#' Filter a vector of taxa to those belonging to a given clade
#'
#' Given a vector of taxon names and a clade name, returns only those taxa
#' whose lineage includes the specified clade.
#'
#' @param taxa A character vector of taxon names.
#' @param clade A character string giving the clade to filter by.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A character vector of taxa that are members of the specified clade.
#'
#' @export
#' @examples
#' \donttest{
#' taxa <- c("Tyrannosaurus", "Triceratops", "Velociraptor",
#'           "Brachiosaurus", "Homo")
#' filter_clade(taxa, "Theropoda")
#' filter_clade(taxa, "Dinosauria")
#' }
filter_clade <- function(taxa, clade, verbose = FALSE) {
  keep <- purrr::map_lgl(taxa, function(t) {
    result <- is_member(t, clade, verbose = verbose)
    if (is.null(result)) FALSE else result
  })
  taxa[keep]
}

#' Plot method for taxodist_cluster objects
#'
#' Plots the hierarchical clustering dendrogram produced by [taxo_cluster()].
#'
#' @param x A `taxodist_cluster` object from [taxo_cluster()].
#' @param main Plot title. Default: `"Taxonomic Clustering"`.
#' @param xlab X-axis label. Default: `""`.
#' @param sub Subtitle. Default: clustering method used.
#' @param ... Additional arguments passed to [stats::plot.hclust()].
#'
#' @return Invisibly returns `x`. Called for side effects (plotting).
#' @export
#' @examples
#' \donttest{
#' taxa <- c("Tyrannosaurus", "Velociraptor", "Triceratops", "Brachiosaurus")
#' cl <- taxo_cluster(taxa)
#' plot(cl)
#' }
plot.taxodist_cluster <- function(x, main = "Taxonomic Clustering",
                                  xlab = "",
                                  sub  = paste("Method:", x$hclust$method),
                                  ...) {
  plot(x$hclust, main = main, xlab = xlab, sub = sub, ...)
  invisible(x)
}

#' Plot method for taxodist_ord objects
#'
#' Plots the PCoA ordination produced by [taxo_ordinate()].
#'
#' @param x A `taxodist_ord` object from [taxo_ordinate()].
#' @param main Plot title. Default: `"Taxonomic Ordination (PCoA)"`.
#' @param xlab X-axis label. Default: `"PC1"`.
#' @param ylab Y-axis label. Default: `"PC2"`.
#' @param labels Character vector of labels. Default: rownames of points.
#' @param ... Additional arguments passed to [graphics::text()].
#'
#' @return Invisibly returns `x`. Called for side effects (plotting).
#' @export
#' @examples
#' \donttest{
#' taxa <- c("Tyrannosaurus", "Velociraptor", "Triceratops", "Brachiosaurus")
#' ord <- taxo_ordinate(taxa)
#' plot(ord)
#' }
plot.taxodist_ord <- function(x, main   = "Taxonomic Ordination (PCoA)",
                              xlab   = "PC1",
                              ylab   = "PC2",
                              labels = rownames(x$points),
                              ...) {
  gof <- round(x$GOF[1], 3)
  graphics::plot(x$points, type = "n",
                 main = paste0(main, "  (GOF = ", gof, ")"),
                 xlab = xlab, ylab = ylab)
  graphics::text(x$points, labels = labels, ...)
  invisible(x)
}

#' Summary method for taxodist_ord objects
#'
#' Prints a summary of the PCoA ordination, including the goodness-of-fit
#' and the proportion of variance explained by each principal coordinate.
#'
#' @param object A \code{taxodist_ord} object from \code{\link{taxo_ordinate}}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns a data frame containing the eigenvalues and
#'   variance explained by each dimension.
#' @export
summary.taxodist_ord <- function(object, ...) {
  cli::cli_h2("Taxonomic Ordination Summary (PCoA)")

  gof <- round(object$GOF[1] * 100, 2)
  cli::cli_text("Goodness-of-Fit (GOF): {gof}%")
  cli::cli_text("")

  if (!is.null(object$eig)) {
    eig_pos <- object$eig[object$eig > 0]
    var_exp <- (eig_pos / sum(eig_pos)) * 100

    k <- ncol(object$points)

    df <- data.frame(
      Axis           = paste0("PC", 1:k),
      Eigenvalue     = object$eig[1:k],
      Variance_Pct   = var_exp[1:k],
      Cumulative_Pct = cumsum(var_exp)[1:k]
    )

    print(df, row.names = FALSE, digits = 4)
    invisible(df)
  } else {
    cli::cli_alert_warning("Eigenvalues not found in the object.")
    invisible(NULL)
  }
}

#' Plot a taxonomic heatmap
#'
#' Computes pairwise taxonomic distances and plots a heatmap with hierarchical
#' clustering dendrograms on the margins. Darker/hotter colors typically
#' represent smaller distances (closer relatives).
#'
#' @param taxa A character vector of taxon names, or a \code{dist} object from
#'   \code{\link{distance_matrix}}.
#' @param ... Additional arguments passed to \code{\link[stats]{heatmap}}.
#'
#' @return Invisibly returns the underlying \code{dist} object.
#'   Called primarily for its side effect (plotting).
#' @export
#' @examples
#' \donttest{
#' taxa <- c("Tyrannosaurus", "Velociraptor", "Homo", "Panthera", "Quercus")
#' taxo_heatmap(taxa)
#' }
taxo_heatmap <- function(taxa, ...) {
  d <- if (inherits(taxa, "dist")) taxa else distance_matrix(taxa, ...)
  mat <- as.matrix(d)
  stats::heatmap(mat, symm = TRUE, margins = c(10, 10), ...)
  invisible(d)
}

#' Get the taxonomic path between two taxa
#'
#' Returns the full node-by-node path from one taxon up to their most recent
#' common ancestor (MRCA) and back down to the other taxon. The result is a
#' data frame with one row per node, making it easy to inspect, filter, or
#' pipe into other functions.
#'
#' @param taxon_a A character string giving the first taxon name.
#' @param taxon_b A character string giving the second taxon name.
#' @param verbose Logical. If `TRUE`, prints progress messages. Default `FALSE`.
#'
#' @return A data frame of class `"taxodist_path"` with columns:
#' \describe{
#'   \item{`node`}{Character. The clade or taxon name at this step.}
#'   \item{`depth`}{Integer. The depth of this node in the full lineage of
#'     its side (or the shared lineage for the MRCA).}
#'   \item{`direction`}{Character. One of `"a"` (ascending from taxon A to
#'     MRCA), `"mrca"` (the shared ancestor), or `"b"` (descending from MRCA
#'     to taxon B).}
#' }
#' Returns `NULL` if either taxon cannot be found.
#'
#' @seealso [mrca()], [shared_clades()], [compare_lineages()]
#' @export
#' @examples
#' \donttest{
#' taxo_path("Tyrannosaurus", "Velociraptor")
#' taxo_path("Tyrannosaurus", "Homo")
#' }
taxo_path <- function(taxon_a, taxon_b, verbose = FALSE) {
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

  result  <- .compute_distance(lin_a, lin_b, taxon_a, taxon_b)
  mrca_d  <- result$mrca_depth

  side_a_nodes <- rev(lin_a[seq_len(mrca_d)])
  side_a <- data.frame(
    node      = side_a_nodes[-length(side_a_nodes)],
    depth     = rev(seq_len(mrca_d - 1L)),
    direction = "a",
    stringsAsFactors = FALSE
  )

  mrca_row <- data.frame(
    node      = result$mrca,
    depth     = mrca_d,
    direction = "mrca",
    stringsAsFactors = FALSE
  )

  side_b_nodes <- lin_b[(mrca_d + 1L):length(lin_b)]
  side_b <- data.frame(
    node      = side_b_nodes,
    depth     = seq_along(side_b_nodes),
    direction = "b",
    stringsAsFactors = FALSE
  )

  path <- rbind(side_a, mrca_row, side_b)
  rownames(path) <- NULL
  structure(path, class = c("taxodist_path", "data.frame"),
            taxon_a = taxon_a, taxon_b = taxon_b)
}

#' Print method for taxodist_path objects
#'
#' @param x A `taxodist_path` object from [taxo_path()].
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns `x`. Called for side effects (printing).
#' @export
print.taxodist_path <- function(x, ...) {
  taxon_a <- attr(x, "taxon_a")
  taxon_b <- attr(x, "taxon_b")
  cli::cli_h2("Taxonomic Path: {taxon_a} \u2192 {taxon_b}")
  for (i in seq_len(nrow(x))) {
    node <- x$node[i]
    dir  <- x$direction[i]
    dep  <- x$depth[i]
    if (dir == "a") {
      cli::cli_bullets(c(" " = "{.emph {node}}  (depth {dep})  \u2191"))
    } else if (dir == "mrca") {
      cli::cli_bullets(c("*" = "{.strong {node}}  (MRCA, depth {dep})"))
    } else {
      cli::cli_bullets(c(" " = "{.emph {node}}  (depth {dep})  \u2193"))
    }
  }
  invisible(x)
}
