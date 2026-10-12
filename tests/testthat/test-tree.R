# ── taxo_tree ────────────────────────────────────────────────────────────────

dinosaur_lineages <- list(
  Tyrannosaurus = c("Biota", "Dinosauria", "Theropoda", "Tyrannosaurus"),
  Velociraptor  = c("Biota", "Dinosauria", "Theropoda", "Velociraptor"),
  Triceratops   = c("Biota", "Dinosauria", "Ornithischia", "Triceratops"),
  Homo          = c("Biota", "Mammalia", "Homo")
)

test_that("taxo_tree builds a phylo object with named clades", {
  tree <- taxo_tree(taxo_from_lineages(dinosaur_lineages), progress = FALSE)
  expect_s3_class(tree, "phylo")
  expect_equal(tree$Nnode, 3L)
  expect_setequal(tree$tip.label, names(dinosaur_lineages))
  expect_equal(tree$node.label, c("Biota", "Dinosauria", "Theropoda"))
  expect_equal(nrow(tree$edge), length(tree$edge.length))
  expect_equal(tree$edge[1, 1], length(dinosaur_lineages) + 1L)
  expect_true(all(tree$edge.length >= 0))
})

test_that("taxo_tree cophenetic distances reproduce distance_matrix", {
  skip_if_not_installed("ape")
  lineages <- c(dinosaur_lineages, list(
    Dinosauria = c("Biota", "Dinosauria"),
    Rex1 = dinosaur_lineages$Tyrannosaurus,
    Rex2 = dinosaur_lineages$Tyrannosaurus,
    Theropoda = c("Biota", "Dinosauria", "Theropoda"),
    Theropoda2 = c("Biota", "Dinosauria", "Theropoda")
  ))
  resolution <- taxo_from_lineages(lineages)
  tree <- taxo_tree(resolution, progress = FALSE)
  expected <- as.matrix(distance_matrix(resolution, progress = FALSE))
  observed <- ape::cophenetic.phylo(tree)[rownames(expected), colnames(expected)]

  expect_equal(observed, expected)
  expect_true(ape::is.ultrametric(tree))
  expect_true(ape::is.rooted(tree))
  expect_no_error(ape::write.tree(tree))
  expect_true("" %in% tree$node.label)
})

test_that("taxo_tree labels a group of identical lineages with their clade", {
  tree <- taxo_tree(taxo_from_lineages(list(
    A = c("Biota", "Alpha"),
    B = c("Biota", "Alpha"),
    C = c("Biota", "Beta")
  )), progress = FALSE)
  expect_equal(tree$node.label, c("Biota", "Alpha"))
  expect_equal(tree$Nnode, 2L)
})

test_that("taxo_tree accepts names and bundles", {
  mockery::stub(taxo_tree, "get_lineage", function(taxon, ...) {
    dinosaur_lineages[[taxon]]
  })
  tree <- taxo_tree(names(dinosaur_lineages), progress = FALSE)
  expect_s3_class(tree, "phylo")

  bundle <- taxo_bundle(taxo_from_lineages(dinosaur_lineages), progress = FALSE)
  expect_equal(taxo_tree(bundle)$Nnode, 3L)
})

test_that("taxo_tree drops taxa without lineages and validates inputs", {
  resolution <- taxo_from_lineages(dinosaur_lineages[1:2])
  resolution$lineage[2] <- list(NULL)
  expect_error(
    expect_warning(taxo_tree(resolution), "without a lineage"),
    "At least two taxa"
  )

  three <- taxo_from_lineages(dinosaur_lineages[1:3])
  three$lineage[3] <- list(NULL)
  expect_warning(tree <- taxo_tree(three), "Triceratops")
  expect_equal(sort(tree$tip.label), c("Tyrannosaurus", "Velociraptor"))

  expect_error(
    taxo_tree(taxo_from_lineages(list(
      A = c("Biota", "A"),
      B = c("Natura", "B")
    ))),
    "same root"
  )
})
