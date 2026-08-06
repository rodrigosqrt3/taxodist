#' Reference Taxonomic Dataset (taxobase)
#'
#' A pre-computed reference dataset for 52 candidate genera distributed across
#' major taxonomic groups. The data were retrieved from The Taxonomicon using
#' the package functions and are intended for examples, documentation, and
#' offline inspection of package output structures.
#'
#' @docType data
#' @name taxobase
#' @usage data(taxobase)
#' @format A structured list containing the following pre-computed components:
#' \describe{
#'   \item{taxa}{Character vector of 52 candidate genera.}
#'   \item{found_taxa}{Character vector of candidates found in The Taxonomicon.}
#'   \item{coverage}{Named logical vector indicating database presence of each taxon.}
#'   \item{matrix}{A pre-computed taxonomic `dist` object for `found_taxa`.}
#'   \item{pairwise}{A pre-computed distance result for *Tyrannosaurus* and *Homo*.}
#'   \item{lineage_homo}{Taxonomic lineage returned for *Homo*.}
#'   \item{lineage_tyrannosaurus}{Taxonomic lineage returned for *Tyrannosaurus*.}
#'   \item{closest}{Result of a closest-relative query for *Tyrannosaurus*.}
#'   \item{filter}{Taxa in `found_taxa` assigned to Theropoda.}
#'   \item{search}{Reference search result for the term `"Bacteria"`.}
#'   \item{statistical_taxa}{Character vector used by the statistical applications vignette.}
#'   \item{statistical_matrix}{Pre-computed `dist` object for `statistical_taxa`.}
#'   \item{metadata}{Data source, generation date, package version, and distance definition.}
#' }
#' @source \url{http://taxonomicon.taxonomy.nl}
#' @examples
#' data(taxobase)
#' length(taxobase$found_taxa)
#' head(labels(taxobase$matrix))
"taxobase"
