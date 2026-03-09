# taxodist

**Taxonomic distance and phylogenetic lineage computation for any taxon on Earth.**

`taxodist` retrieves full hierarchical lineages from [The Taxonomicon](http://taxonomicon.taxonomy.nl) and computes a tree metric distance between any two taxa.

## Installation

```r
# Development version
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

The default metric is the Jaccard-based taxonomic distance:

  $$d_{\text{jaccard}}(A, B) = 1 - \frac{\text{depth}(\text{MRCA}(A,B))}{\text{depth}(A) + \text{depth}(B) - \text{depth}(\text{MRCA}(A,B))}$$

Returns a value between 0 (identical) and 1 (no shared ancestry beyond root),
normalized for lineage depth. Three methods available via `method` argument:
`"jaccard"` (default), `"norm"`, and `"raw"`. The raw metric satisfies the
triangle inequality.

## Data source

All lineage data is sourced from **The Taxonomicon** (taxonomy.nl), based on *Systema Naturae 2000* by S.J. Brands (1989 onwards). The Taxonomicon provides exceptionally deep lineage resolution, substantially exceeding other programmatic sources.

Please cite The Taxonomicon in any published work:

> Brands, S.J. (1989 onwards). *Systema Naturae 2000*. Amsterdam, The Netherlands. Retrieved from The Taxonomicon, http://taxonomicon.taxonomy.nl.
