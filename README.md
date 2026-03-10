# taxodist

[![R-CMD-check](https://github.com/rodrigosqrt3/taxodist/actions/workflows/r.yml/badge.svg)](https://github.com/rodrigosqrt3/taxodist/actions/workflows/r.yml)
[![Coverage](https://codecov.io/gh/rodrigosqrt3/taxodist/branch/main/graph/badge.svg)](https://codecov.io/gh/rodrigosqrt3/taxodist)

**Taxonomic distance and phylogenetic lineage computation for any taxon on Earth.**

`taxodist` retrieves full hierarchical lineages from [The Taxonomicon](http://taxonomicon.taxonomy.nl) and computes a tree metric distance between any two taxa: a pair of dinosaurs, a dinosaur and a fungus, two species of fly, or an oak tree and a human.

## Installation

```r
devtools::install_github("rodrigosqrt3/taxodist")
```

## Basic usage

```r
library(taxodist)

# Get a full lineage
get_lineage("Tyrannosaurus")

# Distance between two taxa
taxo_distance("Tyrannosaurus", "Velociraptor")

# Most recent common ancestor
mrca("Tyrannosaurus", "Triceratops")   # "Dinosauria"
mrca("Tyrannosaurus", "Homo")          # "Amniota"

# Pairwise distance matrix
theropods <- c("Tyrannosaurus", "Velociraptor", "Spinosaurus", "Allosaurus")
distance_matrix(theropods)

# Filter taxa by clade
taxa <- c("Tyrannosaurus", "Triceratops", "Homo", "Quercus")
filter_clade(taxa, "Dinosauria")
```

## The distance metric

`taxodist` measures relatedness by asking a single question: how deep is the most recent common ancestor (MRCA)?

$$d(A, B) = \frac{1}{\text{depth}(\text{MRCA}(A, B))}$$

The deeper the shared ancestor, the smaller the distance and the more related the two taxa are. A shallow MRCA means the two taxa diverged early; a deep MRCA means they share a long common history. The metric returns 0 when one taxon is ancestral to the other, and satisfies the triangle inequality.

The Taxonomicon provides substantially deeper lineage resolution than other programmatic sources, e.g., `Tyrannosaurus` has over 70 nodes in its lineage, which is what makes the distances meaningful across all of life.

## Data source

All lineage data is sourced from **The Taxonomicon** (taxonomy.nl), based on *Systema Naturae 2000* by Sheila J. Brands (1989 onwards). Please cite this resource in any published work using `taxodist`:

> Brands, S.J. (1989 onwards). *Systema Naturae 2000*. Amsterdam, The Netherlands. Retrieved from The Taxonomicon, http://taxonomicon.taxonomy.nl.

## Contributing

Found a taxon with an incorrect lineage? Please [open an issue](https://github.com/rodrigosqrt3/taxodist/issues),
lineage corrections are the most valuable contribution to this package.

