#' Print method for taxodist distance results
#'
#' @param x A list returned by [taxo_distance()].
#' @param ... Additional arguments (ignored).
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
#' \dontrun{
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
#' \dontrun{
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
#' \dontrun{
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
#' \dontrun{
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
