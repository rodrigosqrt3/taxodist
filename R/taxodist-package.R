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
#' The distance metric is based on the depth of the most recent common
#' ancestor (MRCA):
#'
#' \deqn{d(A, B) = \frac{1}{\text{depth}(\text{MRCA}(A,B))}}
#'
#' The deeper the shared ancestor, the smaller the distance. This metric
#' ensures that taxa sharing the same MRCA are always equidistant from any
#' third taxon, regardless of lineage depth below the split — a key
#' biological correctness property absent from Jaccard-based approaches.
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
