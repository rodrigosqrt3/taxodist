# taxodist <picture><source media="(prefers-color-scheme: dark)" srcset="man/figures/taxodist_dark.png"><source media="(prefers-color-scheme: light)" srcset="man/figures/taxodist_sepia.png"><img alt="taxodist logo" src="man/figures/taxodist_sepia.png" align="right" height="200"></picture>

[![CRAN status](https://www.r-pkg.org/badges/version/taxodist)](https://CRAN.R-project.org/package=taxodist) &nbsp; [![R-CMD-check](https://github.com/rodrigosqrt3/taxodist/actions/workflows/r.yml/badge.svg)](https://github.com/rodrigosqrt3/taxodist/actions/workflows/r.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/taxodist/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/taxodist)

**Taxonomic hierarchy distances derived from lineage classifications.**

`taxodist` retrieves taxonomic lineages from [The Taxonomicon](http://taxonomicon.taxonomy.nl) and computes distances from the depth of the most recent common ancestor. It supports comparisons among named taxa at any level represented in the source hierarchy.

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

# Distances between two sets of taxa (rows x columns)
cross_distances(c("Tyrannosaurus", "Triceratops"), c("Velociraptor", "Homo"))

# Classification tree with named clades (an ape-compatible phylo object)
tree <- taxo_tree(theropods)
ape::plot.phylo(tree, show.node.label = TRUE)

# Filter taxa by clade
taxa <- c("Tyrannosaurus", "Triceratops", "Homo", "Quercus")
filter_clade(taxa, "Dinosauria")

# Get the path between two taxa
taxo_path("Tyrannosaurus", "Velociraptor")

# Save and restore the lineage cache across sessions
save_cache()   # persistent per-user location; or save_cache("my_cache.rds")
load_cache()

# Resolve a batch first and preserve ambiguity, coverage, IDs, and lineages;
# `within` restricts homonyms to a clade
resolved <- taxo_resolve(
  c("Tyrannosaurus", "Nereis", "Unknown taxon"),
  within = "Animalia"
)
distance_matrix(resolved)

# Work completely offline with curated or unpublished lineages
local <- taxo_from_lineages(list(
  Alpha = c("Biota", "Animalia", "Alpha"),
  Beta = c("Biota", "Animalia", "Beta")
), source = "My curated taxonomy")
distance_matrix(local)

# Preserve the complete analysis and exchange it across implementations
bundle <- taxo_bundle(resolved)
write_taxodist_bundle(bundle, "analysis.json")
restored <- read_taxodist_bundle("analysis.json")

# Months later: what changed in the source classification?
compare_bundles("analysis.json", taxo_bundle(restored$resolution$input))
```

## The distance metric

`taxodist` measures relatedness by asking a single question: how deep is the most recent common ancestor (MRCA)?

$$
d(A,B) =
\begin{cases}
0, & A = B, \\
\dfrac{1}{\text{depth}(\text{MRCA}(A,B))}, & A \ne B.
\end{cases}
$$

The deeper the shared ancestor, the smaller the distance and the more related the two taxa are. A shallow MRCA means the two taxa diverged early; a deep MRCA means they share a long common history. Zero is reserved for identical hierarchy nodes. Consequently, a taxon and one of its descendants have a positive distance even though they are connected by ancestry. This distinction makes the measure a proper ultrametric on each connected hierarchy.

Distance and membership answer different questions. For example, *Tyrannosaurus* has a positive distance from *Dinosauria* because they are distinct nodes, while `is_member("Tyrannosaurus", "Dinosauria")` returns `TRUE`. Use `is_member()` or `taxo_path()` when the relationship of interest is containment or ancestry.

The numerical values depend on the resolution of the classification returned by The Taxonomicon. They represent separation within that hierarchy and should not be interpreted as evolutionary time or phylogenetic branch length.

## Caching

Lineages are cached in memory automatically during a session. To persist the
cache across sessions and avoid redundant network requests, use `save_cache()`
and `load_cache()`, which default to a per-user cache directory, or pass an
explicit file path. Requests are retried automatically after transient network
errors (`options(taxodist.max_tries = ...)`, default 3).

## Data source

All lineage data is sourced from **The Taxonomicon** (taxonomy.nl), based on *Systema Naturae 2000* by Sheila J. Brands (1989 onwards). Please cite this resource in any published work using `taxodist`:

> Brands, S.J. (1989 onwards). *Systema Naturae 2000*. Amsterdam, The Netherlands. Retrieved from The Taxonomicon, http://taxonomicon.taxonomy.nl.

## Contributing

Found a taxon with an incorrect lineage? Please [open an issue](https://github.com/rodrigosqrt3/taxodist/issues),
lineage corrections are the most valuable contribution to this package.
