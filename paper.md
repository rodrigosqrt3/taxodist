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
`taxodist` retrieves their full lineages, identifies their most recent common 
ancestor (MRCA), and returns a dissimilarity index reflecting their shared 
evolutionary history. The package supports individual distance queries, pairwise 
distance matrices, clade membership testing, lineage comparison utilities, and 
session-level caching to minimise redundant network requests.

# Statement of Need

Quantifying relatedness between taxa is a fundamental operation in comparative
biology, ecology, and paleontology. Existing programmatic solutions such as
`ape` [@paradis2004], `taxize` [@chamberlain2013], and the Open Tree of Life
[@hinchliff2015] provide access to taxonomic data and phylogenetic tools, but
computing a pairwise dissimilarity between arbitrary taxa, without requiring
a pre-built phylogenetic tree object, remains cumbersome. A researcher wishing
to quantify how related *Tyrannosaurus* is to *Homo sapiens*, or to rank a list
of candidate taxa by relatedness to a query taxon, must either construct a
formal phylogeny or rely on manual inspection of lineage data.

`taxodist` fills this gap by providing a lightweight, tree-free interface for
taxonomic distance computation. It requires only taxon names as input and
returns a dissimilarity index that correctly reflects the depth of shared
ancestry. The package is particularly useful for exploratory analyses,
educational contexts, and research in paleontology and macroecology where
formal phylogenetic trees may not be available or necessary.

# State of the Field

Several R packages address related needs. `ape` [@paradis2004] is the
standard tool for phylogenetic analysis in R, but requires a formal `phylo`
object constructed from sequence data or a pre-existing tree. `taxize`
[@chamberlain2013] provides programmatic access to multiple taxonomic
databases and can retrieve lineage information, but does not compute pairwise
dissimilarity indices directly. 

`taxodist` differs from these tools in two key ways. First, it requires no
pre-built tree: distances are computed directly from hierarchical lineage
strings retrieved at query time. Second, it uses The Taxonomicon as its data
source, which provides substantially deeper lineage resolution than other
programmatic sources. For example, the genus *Tyrannosaurus* has over 70
nodes in its lineage in The Taxonomicon.

# Software Design

The core design decision in `taxodist` was to define a dissimilarity index
based solely on the depth of the MRCA, rather than on the full path length
between taxa:

$$d(A, B) = \frac{1}{\text{depth}(\text{MRCA}(A, B))}$$

This choice was made deliberately to achieve independence from post-divergence
lineage resolution. A path-length metric, counting edges between two taxa on
the tree, penalises taxa whose lineages are more finely resolved, producing
artifactually larger distances for well-studied groups. The MRCA-depth index
avoids this: two taxa that diverged at the same node are always equidistant
from any third taxon, regardless of how many intermediate clades have been
described below the split point.

The tradeoff is that the index is a dissimilarity measure rather than a
formal metric in the mathematical sense, as the triangle inequality is not
guaranteed for all possible lineage configurations. This is documented
explicitly in the package.

A session-level in-memory cache is used to avoid redundant network requests.
Lineage retrieval involves HTML parsing of The Taxonomicon's tree pages, which
required careful handling of navigation elements, *incertae sedis* child taxa,
and rank-prefix strings that appear alongside the true ancestor chain. The
parser truncates the link list at the first child-level entry, ensuring that
only genuine ancestral nodes are included in the returned lineage.

# Research Impact Statement

`taxodist` was developed to address a concrete gap identified during research
in vertebrate paleontology and macroecology, where pairwise taxonomic
distances between taxa from diverse groups (including extinct and extant
clades) are needed without the overhead of constructing formal phylogenies.
The package is currently in its initial release (version 0.1.0) and has not
yet been cited in published work. However, it addresses a gap that is not
filled by any existing R package, and its coverage of all of Biota, makes it 
applicable across a wide range of biological disciplines. The package passes 
R CMD check on Windows, macOS, and Ubuntu across four R versions, has 100% test 
coverage, and is documented with a full vignette.

# AI Usage Disclosure

Generative AI tools were used during the development of `taxodist` for iterative debugging 
of the HTML parser, test suite development. All AI-generated code and text was 
reviewed, tested, and verified by the author. The dissimilarity index, the scraping 
architecture, and the scientific design decisions are the author's own. No AI-generated 
content was used without verification against primary sources or direct empirical testing.

# Acknowledgements

The author thanks Sheila J. Brands for maintaining The Taxonomicon, an
extraordinary curated resource that makes this work possible, and for her
decades of dedicated work on *Systema Naturae 2000*. This work received no
financial support.

# References


