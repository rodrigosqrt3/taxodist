# taxodist <picture><source media="(prefers-color-scheme: dark)" srcset="man/figures/taxodist_dark.png"><source media="(prefers-color-scheme: light)" srcset="man/figures/taxodist_sepia.png"><img alt="taxodist logo" src="man/figures/taxodist_sepia.png" align="right" height="200"></picture>

[![CRAN status](https://www.r-pkg.org/badges/version/taxodist)](https://CRAN.R-project.org/package=taxodist) &nbsp; [![R-CMD-check](https://github.com/rodrigosqrt3/taxodist/actions/workflows/r.yml/badge.svg)](https://github.com/rodrigosqrt3/taxodist/actions/workflows/r.yml) &nbsp; [![codecov](https://codecov.io/gh/rodrigosqrt3/taxodist/branch/main/graph/badge.svg)](https://app.codecov.io/gh/rodrigosqrt3/taxodist)

**Taxonomic hierarchy distances derived from lineage classifications.**

`taxodist` retrieves taxonomic lineages and computes distances from the depth of the most recent common ancestor. Lineages come from [The Taxonomicon](http://taxonomicon.taxonomy.nl) by default, a classification with very deep lineages and extensive fossil coverage, or from the [Open Tree of Life](https://tree.opentreeoflife.org/about/taxonomy-version) reference taxonomy or the [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy). It supports comparisons among named taxa at any level represented in the source hierarchy, offers alternative distance metrics, measures how much results depend on the classification, and provides community analyses: Clarke and Warwick's taxonomic distinctness, taxonomic dissimilarity between samples, and standardized effect sizes (NRI, NTI).

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

## Taxonomic sources

```r
taxo_sources()   # sources, identifiers, licences, and citations

# The same taxa in the Open Tree Taxonomy or the NCBI Taxonomy
apes <- c("Homo", "Pan", "Gorilla", "Pongo")
ott <- taxo_resolve(apes, source = "ott")     # identifiers such as "ott:770309"
ncbi <- taxo_resolve(apes, source = "ncbi")   # identifiers such as "ncbi:9605"
distance_matrix(ott)                          # every matrix, tree, and bundle function works

# How much do the distances depend on the classification?
compare_sources(c("Homo", "Pan", "Mus", "Gallus", "Danio", "Quercus"),
                sources = c("taxonomicon", "ott", "ncbi"))
```

Distances are only comparable within one source: The Taxonomicon resolves many more unranked clades, so its lineages are typically two to three times deeper. Resolutions and bundles record the source and, for the Open Tree Taxonomy, its version. Clade names for `within` follow each source (animals are `"Animalia"` in The Taxonomicon and `"Metazoa"` in the others).

## Community analyses

```r
# comm: samples in rows, taxa in columns (presence or abundance)
d <- distance_matrix(colnames(comm), metric = "ranked")

taxo_distinctness(comm, d)                    # Delta, Delta*, Delta+, Lambda+
taxo_dissimilarity(comm, d)                   # Gamma+ between samples, a dist object
taxo_ses(comm, d, index = "mpd")              # NRI against a species-pool null model
taxo_ses(comm, d, index = "mntd", null_model = "independent_swap")   # NTI
```

`taxo_dissimilarity()` compares samples even when they share no taxon (Clarke, Somerfield & Chapman 2006) and reduces to the Sørensen and Kulczynski dissimilarities on a flat classification. `taxo_ses()` follows Webb et al. (2002) and the definitions of the picante package. See `vignette("sources-and-communities")`.

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

The numerical values depend on the resolution of the source classification. They represent separation within that hierarchy and should not be interpreted as evolutionary time or phylogenetic branch length.

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
requests to each source are spaced by at least 0.25 s (`options(taxodist.min_interval = ...)`)
so that large batches do not overload the servers. NCBI requests are spaced by
at least 0.4 s, below NCBI's limit of three per second; with an
[NCBI API key](https://www.ncbi.nlm.nih.gov/books/NBK25497/) in the
`NCBI_API_KEY` environment variable the limit rises to ten per second.
`taxonomicon_available()` checks quickly whether The Taxonomicon can be reached.

## Data sources

Please cite the source of the lineages in any published work using `taxodist` (`citation("taxodist")` gives all entries):

- **The Taxonomicon** (default): Brands, S.J. (1989 onwards). *Systema Naturae 2000*. Universal Taxonomic Services, Zwaag, The Netherlands. http://taxonomicon.taxonomy.nl.
- **Open Tree Taxonomy** (`source = "ott"`): Rees, J.A. & Cranston, K. (2017). Automated assembly of a reference taxonomy for phylogenetic data synthesis. *Biodiversity Data Journal* 5: e12581. doi:10.3897/BDJ.5.e12581.
- **NCBI Taxonomy** (`source = "ncbi"`): Schoch, C.L. et al. (2020). NCBI Taxonomy: a comprehensive update on curation, resources and tools. *Database* 2020: baaa062. doi:10.1093/database/baaa062.

## Contributing

Found a taxon with an incorrect lineage? Please [open an issue](https://github.com/rodrigosqrt3/taxodist/issues),
lineage corrections are the most valuable contribution to this package.
