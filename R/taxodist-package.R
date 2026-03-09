#' taxodist: Taxonomic Distance and Phylogenetic Lineage Computation
#'
#' @description
#' taxodist computes phylogenetic distances between any two taxa using
#' hierarchical lineage data retrieved from The Taxonomicon
#' (taxonomy.nl), a comprehensive curated classification of all life
#' based on Systema Naturae 2000.
#'
#' ## Core functions
#'
#' - [get_lineage()] — retrieve the full lineage of any taxon
#' - [taxo_distance()] — compute the tree metric distance between two taxa
#' - [mrca()] — find the most recent common ancestor
#' - [distance_matrix()] — compute all pairwise distances for a set of taxa
#' - [closest_relative()] — find the closest relative among candidates
#' - [compare_lineages()] — print a side-by-side lineage comparison
#' - [shared_clades()] — list clades shared between two taxa
#' - [is_member()] — test clade membership
#' - [filter_clade()] — filter taxa by clade membership
#' - [check_coverage()] — check Taxonomicon coverage for a list of taxa
#' - [lineage_depth()] — get the lineage depth of a taxon
#' - [clear_cache()] — clear the session lineage cache
#'
#' ## Mathematical background
#'
#' The default metric is the Jaccard-based taxonomic distance:
#'
#' \deqn{d_{jaccard}(A, B) = 1 - \frac{\text{depth}(\text{MRCA}(A,B))}
#' {\text{depth}(A) + \text{depth}(B) - \text{depth}(\text{MRCA}(A,B))}}
#'
#' This returns a value in \eqn{[0,1]}, normalized for lineage depth so that
#' deeply resolved taxa are not penalized. The underlying raw metric
#' \eqn{d_{raw}(A,B) = \text{depth}(A) + \text{depth}(B) -
#' 2 \cdot \text{depth}(\text{MRCA}(A,B))} satisfies the triangle inequality.
#' Three methods are available: \code{"jaccard"} (default), \code{"norm"},
#' and \code{"raw"}.
#'
#' ## Data source
#'
#' All lineage data is sourced from The Taxonomicon (taxonomy.nl), based on
#' Systema Naturae 2000 by S.J. Brands (1989 onwards). Please cite this
#' resource when using taxodist in published work.
#'
#' @references
#' Brands, S.J. (1989 onwards). Systema Naturae 2000. Amsterdam,
#' The Netherlands. Retrieved from The Taxonomicon,
#' \url{http://taxonomicon.taxonomy.nl}.
#'
#' @docType package
#' @name taxodist-package
#' @aliases taxodist
"_PACKAGE"
