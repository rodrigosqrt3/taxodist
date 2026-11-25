#' taxodist: Taxonomic Hierarchy Distances and Lineage Analysis
#'
#' @description
#' taxodist retrieves taxonomic lineages from The Taxonomicon (taxonomy.nl),
#' the Open Tree of Life reference taxonomy, or the NCBI Taxonomy, or accepts
#' user-supplied lineages, then computes distances from the depth of the most
#' recent common ancestor.
#' The resulting values describe separation within the source classification;
#' they are not estimates of evolutionary time or phylogenetic branch length.
#'
#' ## Core functions
#'
#' - [get_lineage()] — retrieve the full lineage of any taxon
#' - [taxo_resolve()] — resolve and audit a batch of names or IDs, in any
#'   source
#' - [taxo_sources()] — list the supported taxonomic sources
#' - [compare_sources()] — compare the distances given by several sources
#' - [taxo_from_lineages()] — create an offline resolution from local lineages
#' - [taxo_distance()] — compute the ultrametric distance between two taxa
#' - [mrca()] — find the most recent common ancestor
#' - [distance_matrix()] — compute all pairwise distances for a set of taxa
#' - [cross_distances()] — compute distances between two sets of taxa
#' - [taxo_tree()] — build a `phylo` tree whose nodes are named clades
#' - [closest_relative()] — find the closest relative among candidates
#' - [compare_lineages()] — print a side-by-side lineage comparison
#' - [shared_clades()] — list clades shared between two taxa
#' - [is_member()] — test clade membership
#' - [filter_clade()] — filter taxa by clade membership
#' - [check_coverage()] — check Taxonomicon coverage for a list of taxa
#' - [lineage_depth()] — get the lineage depth of a taxon
#' - [taxo_distinctness()] / [taxo_distinctness_test()] — taxonomic
#'   distinctness of samples
#' - [taxo_dissimilarity()] — taxonomic dissimilarity between samples
#' - [taxo_ses()] — standardized effect sizes of community structure (NRI, NTI)
#' - [taxo_bundle()] — preserve an analysis and its provenance
#' - [compare_bundles()] — report what changed between two analyses
#' - [save_cache()] / [load_cache()] — persist retrieved lineages
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
#' ## Data sources
#'
#' Online lineage data is sourced from The Taxonomicon (taxonomy.nl), based on
#' Systema Naturae 2000 by S.J. Brands (1989 onwards), by default; from the
#' Open Tree Taxonomy (Rees and Cranston 2017) with `source = "ott"`; or from
#' the NCBI Taxonomy (Schoch et al. 2020) with `source = "ncbi"`. Please cite
#' the source used (see [taxo_sources()] and `citation("taxodist")`), or record
#' the relevant source when using [taxo_from_lineages()].
#'
#' ## Options
#'
#' - `taxodist.max_tries`: number of attempts made for each request (default
#'   `3`). Network errors, server errors (5xx), and refused requests (429) are
#'   retried with exponential backoff; other client errors (4xx) are not.
#' - `taxodist.min_interval`: minimum time in seconds between consecutive
#'   requests to each source (default `0.25`), so that large batches do not
#'   overload the servers. Set to `0` to disable. NCBI requests are always
#'   spaced by at least 0.4 s (0.12 s with an API key).
#' - `taxodist.ncbi_api_key`: NCBI API key; defaults to the `NCBI_API_KEY`
#'   environment variable.
#' - `taxodist.rank_levels`: ranks used by the `"ranked"` metric.
#'
#' @references
#' Brands, S.J. (1989 onwards). Systema Naturae 2000. Amsterdam,
#' The Netherlands. Retrieved from The Taxonomicon,
#' \url{http://taxonomicon.taxonomy.nl}.
#'
#' Rees, J.A. and Cranston, K. (2017). Automated assembly of a reference
#' taxonomy for phylogenetic data synthesis. *Biodiversity Data Journal*, 5,
#' e12581. \doi{10.3897/BDJ.5.e12581}
#'
#' Schoch, C.L. et al. (2020). NCBI Taxonomy: a comprehensive update on
#' curation, resources and tools. *Database*, 2020, baaa062.
#' \doi{10.1093/database/baaa062}
#'
#' @name taxodist-package
#' @aliases taxodist
"_PACKAGE"
