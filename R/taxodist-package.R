#' taxodist: Taxonomic Hierarchy Distances and Lineage Analysis
#'
#' @description
#' taxodist retrieves taxonomic lineages from The Taxonomicon (taxonomy.nl)
#' and computes distances from the depth of the most recent common ancestor.
#' The resulting values describe separation within the source classification;
#' they are not estimates of evolutionary time or phylogenetic branch length.
#'
#' ## Core functions
#'
#' - [get_lineage()] — retrieve the full lineage of any taxon
#' - [taxo_distance()] — compute the ultrametric distance between two taxa
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
#' \deqn{d(A,A) = 0}{d(A,A) = 0}
#' \deqn{d(A,B) = 1 / h(A,B)}{d(A,B) = 1 / h(A,B)}
#'
#' Here, `h(A,B)` denotes the depth of the MRCA. The deeper the shared
#' ancestor, the smaller the distance. Zero is reserved
#' for identical hierarchy nodes; distinct ancestor-descendant pairs therefore
#' have positive distance. The measure is an ultrametric on each connected
#' taxonomic hierarchy. Its numerical values depend on the resolution of the
#' classification returned by the data source.
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
