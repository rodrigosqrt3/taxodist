---
title: 'taxodist: Taxonomic Distance and Phylogenetic Lineage Computation in R'
tags:
  - R
  - taxonomy
  - phylogenetics
  - biodiversity
  - paleontology
authors:
  - name: Rodrigo Fonseca Villa
    orcid: 0009-0005-2938-2270
    affiliation: 1
affiliations:
  - name: Universidade Federal do Rio Grande do Sul (UFRGS), Brazil
    index: 1
date: 11 March 2026
bibliography: paper.bib
---

# Summary

`taxodist` is an R package that computes dissimilarity indices between any two
taxa using hierarchical lineage data retrieved from The Taxonomicon
[@brands1989], a comprehensive curated classification of all life based on
*Systema Naturae 2000*. Given any two taxon names, at any taxonomic level,
from species to order, `taxodist` retrieves their full lineages, identifies
their most recent common ancestor (MRCA), and returns a dissimilarity index
reflecting their shared evolutionary history. The package supports individual
distance queries, pairwise distance matrices, clade membership testing,
lineage comparison utilities, and session-level caching to minimise redundant
network requests.

# Statement of Need

Quantifying relatedness between taxa is a fundamental operation in comparative
biology, ecology, and paleontology. Existing programmatic solutions such as
`ape` [@paradis2004], `taxize` [@chamberlain2013], and the Open Tree of Life
[@hinchliff2015] provide access to taxonomic data and phylogenetic tools, but
computing a pairwise dissimilarity between arbitrary taxa — without requiring
a pre-built phylogenetic tree object — remains cumbersome. A researcher wishing
to quantify how related *Tyrannosaurus* is to *Homo sapiens*, or to rank a list
of candidate taxa by relatedness to a query taxon, must either construct a
formal phylogeny or rely on manual inspection of lineage data.

`taxodist` fills this gap by providing a lightweight, tree-free interface for
taxonomic distance computation. It requires only taxon names as input and
returns a dissimilarity index that correctly reflects the depth of shared
ancestry. Crucially, the index used by `taxodist` is independent of
post-divergence lineage resolution: two taxa that diverged at the same node are
always equidistant from any third taxon, regardless of how many intermediate
clades have been described below the split point. This property makes
`taxodist` particularly well-suited for comparisons across groups with uneven
taxonomic resolution, such as fossil and extant clades.

The Taxonomicon provides substantially deeper lineage resolution than other
programmatic sources such as the Open Tree of Life [@hinchliff2015], for example, the genus *Tyrannosaurus* has over 70 nodes
in its lineage in The Taxonomicon, compared to far fewer in sources such as
NCBI Taxonomy or the Open Tree of Life.

# Methodology

For two taxa A and B with lineages $L_A$ and $L_B$ ordered from root to tip,
`taxodist` identifies the most recent common ancestor as the deepest node
shared by both lineages:

$$\text{MRCA}(A, B) = \underset{n \,\in\, L_A \cap L_B}{\arg\max}\; \text{depth}(n)$$

The dissimilarity index is then defined as:

$$d(A, B) = \frac{1}{\text{depth}(\text{MRCA}(A, B))}$$

where depth is measured as the position of the MRCA node in the lineage ordered
from root to tip. The index returns 0 when one taxon is directly ancestral to
the other, and decreases toward 0 as the MRCA becomes deeper (more recent
in the classification hierarchy). When two taxa share only the root node Biota
at depth 1, the index returns its maximum value of 1. The index is symmetric
and respects the intuition that more deeply nested shared ancestry implies
greater relatedness.

Lineages are retrieved by querying The Taxonomicon's web interface and parsing
the resulting HTML. The parser filters navigation elements and *incertae sedis*
child taxa that appear in the page alongside the true ancestor chain, ensuring
that only genuine ancestral nodes are included in the lineage. A session-level
cache prevents redundant network requests when the same taxon is queried
multiple times within an R session.

# Acknowledgements

The author thanks Sheila J. Brands for maintaining The Taxonomicon, an
extraordinary curated resource that makes this work possible, and for her
decades of dedicated work on *Systema Naturae 2000*.

# References
