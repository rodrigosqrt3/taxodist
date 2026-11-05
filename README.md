# taxodist <picture><source media="(prefers-color-scheme: dark)" srcset="man/figures/taxodist_dark.png"><source media="(prefers-color-scheme: light)" srcset="man/figures/taxodist_sepia.png"><img alt="taxodist logo" src="man/figures/taxodist_sepia.png" align="right" height="200"></picture>

[![CRAN status](https://www.r-pkg.org/badges/version/taxodist)](https://CRAN.R-project.org/package=taxodist) &nbsp; [![R-CMD-check](https://github.com/rodrigosqrt3/taxodist/actions/workflows/r.yml/badge.svg)](https://github.com/rodrigosqrt3/taxodist/actions/workflows/r.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/taxodist/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/taxodist)

**Taxonomic hierarchy distances derived from lineage classifications.**

`taxodist` retrieves taxonomic lineages from [The Taxonomicon](http://taxonomicon.taxonomy.nl) and computes distances from the depth of the most recent common ancestor. It supports comparisons among named taxa at any level represented in the source hierarchy, offers alternative distance metrics, and computes Clarke and Warwick's taxonomic distinctness indices.

The same package is available for [Python](https://github.com/rodrigosqrt3/taxodist-py) and [Julia](https://github.com/rodrigosqrt3/Taxodist.jl); the three implementations give identical results and exchange analyses through a common JSON bundle format.

## Installation

```r
# Release version from CRAN
install.packages("taxodist")

# Development version from GitHub
devtools::install_github("rodrigosqrt3/taxodist")
```

## Basic usage

```r
library(taxodist)

# Get a full lineage, and the rank of each of its nodes
get_lineage("Tyrannosaurus")
lineage_ranks("Tyrannosaurus")

# Distance between two taxa
taxo_distance("Tyrannosaurus", "Velociraptor")

# The same distance with another metric (see "Distance metrics" below)
taxo_distance("Tyrannosaurus", "Velociraptor", metric = "wu_palmer")

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

# Taxonomic distinctness of a community (Clarke & Warwick), and a test
# against random samples from a master list
pool <- c("Tyrannosaurus", "Velociraptor", "Triceratops", "Stegosaurus",
          "Brachiosaurus", "Diplodocus", "Allosaurus", "Spinosaurus")
d <- distance_matrix(pool, metric = "ranked")
taxo_distinctness(c("Tyrannosaurus", "Velociraptor", "Allosaurus"), d)
taxo_distinctness_test(c("Tyrannosaurus", "Velociraptor", "Allosaurus"), d)

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

## Distance metrics

By default, `taxodist` measures relatedness by asking a single question: how deep is the most recent common ancestor (MRCA)?

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

Every distance function accepts a `metric` argument. With *h* the depth of the MRCA and *d<sub>A</sub>*, *d<sub>B</sub>* the depths of the two taxa:

| `metric` | Distance | Properties |
|---|---|---|
| `"reciprocal"` (default) | 1 / *h* | ultrametric; best tracks divergence time within a clade |
| `"path"` | *d<sub>A</sub>* + *d<sub>B</sub>* − 2*h* | number of edges between the two taxa |
| `"wu_palmer"` | 1 − 2*h* / (*d<sub>A</sub>* + *d<sub>B</sub>*) | between 0 and 1; corrects for uneven resolution of the classification |
| `"ranked"` | half the path through the standard ranks only | the step weights of Clarke and Warwick's taxonomic distinctness |

All metrics are 0 for identical taxa. The ranks used by `"ranked"` (kingdom to species by default) can be changed with `options(taxodist.rank_levels = ...)`.

## Caching

Lineages are cached in memory automatically during a session. To persist the
cache across sessions and avoid redundant network requests, use `save_cache()`
and `load_cache()`, which default to a per-user cache directory, or pass an
explicit file path. Requests are retried automatically after transient network
errors (`options(taxodist.max_tries = ...)`, default 3), and consecutive
requests are spaced by at least 0.25 s (`options(taxodist.min_interval = ...)`)
so that large batches do not overload the server. `taxonomicon_available()`
checks quickly whether the server can be reached.

## Data source

All lineage data is sourced from **The Taxonomicon** (taxonomy.nl), based on *Systema Naturae 2000* by Sheila J. Brands (1989 onwards). Please cite this resource in any published work using `taxodist`:

> Brands, S.J. (1989 onwards). *Systema Naturae 2000*. Amsterdam, The Netherlands. Retrieved from The Taxonomicon, http://taxonomicon.taxonomy.nl.

## Contributing

Found a taxon with an incorrect lineage? Please [open an issue](https://github.com/rodrigosqrt3/taxodist/issues),
lineage corrections are the most valuable contribution to this package.
