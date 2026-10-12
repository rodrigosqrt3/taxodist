#' Build a taxonomic tree from lineages
#'
#' Builds a tree whose internal nodes are the named clades of the source
#' classification, directly from the lineages of the supplied taxa. Unlike
#' [taxo_cluster()], which reconstructs a dendrogram from distances, the
#' topology here is the classification itself and internal nodes carry clade
#' names.
#'
#' @param taxa A character vector of taxon names, a `taxodist_resolution`
#'   from [taxo_resolve()] or [taxo_from_lineages()], or a `taxodist_bundle`.
#' @param verbose,progress Passed to [get_lineage()] retrieval when `taxa`
#'   is a character vector.
#'
#' @return An object of class `"phylo"`, compatible with the \pkg{ape}
#'   package (for example `ape::plot.phylo()`, `ape::write.tree()`, and
#'   `ape::cophenetic.phylo()`). Tip labels are the input labels and node
#'   labels are clade names. Taxa without a lineage are dropped with a
#'   warning.
#'
#' @details
#' Clades with a single descendant among the supplied taxa are collapsed, so
#' every internal node is a point where the taxa diverge and is labeled with
#' the deepest clade shared by its descendants (their MRCA).
#'
#' Branch lengths are chosen so that the tree is ultrametric and its
#' cophenetic (patristic) distances reproduce [distance_matrix()] exactly:
#' an internal node for a clade at depth `h` is placed at height `1 / (2h)`
#' and tips at height zero. A taxon that is an ancestor of another supplied
#' taxon is attached to its clade node by a branch of length `1 / (2h)`, and
#' taxa with identical lineages are joined by zero-length branches.
#'
#' All lineages must share the same root node; otherwise distances between
#' some taxa are infinite and no single tree exists.
#'
#' @seealso [taxo_cluster()], [distance_matrix()]
#' @export
#' @examples
#' resolved <- taxo_from_lineages(list(
#'   Tyrannosaurus = c("Biota", "Dinosauria", "Theropoda", "Tyrannosaurus"),
#'   Velociraptor  = c("Biota", "Dinosauria", "Theropoda", "Velociraptor"),
#'   Triceratops   = c("Biota", "Dinosauria", "Ornithischia", "Triceratops"),
#'   Homo          = c("Biota", "Mammalia", "Homo")
#' ))
#' tree <- taxo_tree(resolved)
#' tree$node.label
#' if (requireNamespace("ape", quietly = TRUE)) {
#'   plot(tree, show.node.label = TRUE)
#'   ape::write.tree(tree)
#' }
taxo_tree <- function(taxa, verbose = FALSE, progress = TRUE) {
  input <- .matrix_input(taxa, get_lineage, verbose = verbose, progress = progress)
  labels <- input$labels
  lineages <- input$lineages

  missing <- vapply(lineages, is.null, logical(1))
  if (any(missing)) {
    cli::cli_warn(c(
      "!" = "{sum(missing)} tax{?on/a} without a lineage {?was/were} dropped from the tree.",
      "i" = "{.val {labels[missing]}}"
    ))
    labels <- labels[!missing]
    lineages <- lineages[!missing]
  }
  if (length(lineages) < 2L) {
    cli::cli_abort("At least two taxa with lineages are required to build a tree.")
  }
  lens <- lengths(lineages)
  if (any(lens == 0L) ||
      length(unique(vapply(lineages, `[[`, character(1), 1L))) > 1L) {
    cli::cli_abort(c(
      "x" = "All lineages must share the same root node.",
      "i" = "Taxa without a common ancestor have infinite distance and cannot be placed in one tree."
    ))
  }

  # Recursively group taxa sharing the lineage prefix of length `depth`.
  # Returns either list(tip = <index>) or an internal node
  # list(label, height, children).
  build <- function(idx, depth) {
    clade <- lineages[[idx[[1L]]]][[depth]]
    ended <- idx[lens[idx] == depth]
    continuing <- idx[lens[idx] > depth]

    children <- list()
    if (length(ended) == 1L) {
      children <- list(list(tip = ended))
    } else if (length(ended) > 1L) {
      # Identical lineages: join at height zero so their distance is 0.
      children <- list(list(
        label = if (length(continuing) > 0L) "" else clade,
        height = 0,
        children = lapply(ended, function(i) list(tip = i))
      ))
    }
    if (length(continuing) > 0L) {
      next_nodes <- vapply(
        lineages[continuing], `[[`, character(1), depth + 1L
      )
      groups <- split(continuing, factor(next_nodes, levels = unique(next_nodes)))
      children <- c(children, unname(lapply(groups, build, depth = depth + 1L)))
    }

    if (length(children) == 1L) {
      return(children[[1L]])
    }
    list(label = clade, height = 1 / (2 * depth), children = children)
  }
  root <- build(seq_along(lineages), 1L)

  n_tips <- length(lineages)
  tip_label <- character(n_tips)
  node_label <- character(0)
  edge_parent <- integer(0)
  edge_child <- integer(0)
  edge_length <- numeric(0)
  next_tip <- 0L
  next_node <- n_tips

  # Preorder traversal numbering tips 1..n and internal nodes from n + 1,
  # which yields an edge matrix in "cladewise" order.
  visit <- function(node, parent = NULL, parent_height = NULL) {
    if (!is.null(node$tip)) {
      next_tip <<- next_tip + 1L
      id <- next_tip
      tip_label[[id]] <<- labels[[node$tip]]
      height <- 0
    } else {
      next_node <<- next_node + 1L
      id <- next_node
      node_label[[id - n_tips]] <<- node$label
      height <- node$height
    }
    if (!is.null(parent)) {
      edge_parent <<- c(edge_parent, parent)
      edge_child <<- c(edge_child, id)
      edge_length <<- c(edge_length, parent_height - height)
    }
    for (child in node$children) {
      visit(child, id, height)
    }
  }
  visit(root)

  structure(
    list(
      edge = cbind(edge_parent, edge_child, deparse.level = 0),
      edge.length = edge_length,
      tip.label = tip_label,
      Nnode = next_node - n_tips,
      node.label = node_label
    ),
    class = "phylo",
    order = "cladewise"
  )
}
