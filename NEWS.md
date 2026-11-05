# taxodist 0.1.0

* Initial release
* Compute taxonomic distances between any two taxa using The Taxonomicon
* Session-level caching for lineage data to minimize network requests
* Functions: `taxo_distance()`, `mrca()`, `distance_matrix()`, `closest_relative()`,
  `compare_lineages()`, `shared_clades()`, `is_member()`, `filter_clade()`,
  `lineage_depth()`, `check_coverage()`, `clear_cache()`

# taxodist 0.2.0

## New functions

* `taxo_cluster()`: hierarchical clustering of taxa by taxonomic distance.
* `taxo_ordinate()`: ordination (PCoA) of taxa in taxonomic distance space.
* `save_cache()`: serialises the session lineage cache to an `.rds` file.
* `load_cache()`: restores a previously saved cache, avoiding redundant network requests.
* `taxo_path()`: returns the full node-by-node path between two taxa as a
  tidy data frame, ascending from taxon A to their MRCA and descending to
  taxon B.

## Documentation

* Added a new vignette (`Statistical Applications of taxodist`) demonstrating 
  the integration of `taxodist` with `ape` (tree plotting) and `vegan` 
  (taxonomic distinctness, Mantel tests, and PERMANOVA).

## Minor improvements

* Added `vegan` to `Suggests` in `DESCRIPTION`.
* Documented the compatibility of `distance_matrix()` output with `vegan` 
  functions (`taxondive()`, `mantel()`, `adonis2()`).
* Documented the conversion of `taxo_cluster()` results into `phylo` objects 
  using `ape::as.phylo()`.

## Bug fixes

* Fixed incorrect MRCA computation caused by unnamed crown nodes (`[crown]`)
  being matched across different lineages as if they were the same ancestor.
* Fixed parsing of Taxonomicon rank prefixes (`Cohort`, `Subcohort`, 
  `Magnorder`, `Grandorder`, `Parvorder`, `Legion`) that were being retained 
  as bare rank names instead of the actual clade names (e.g., `Placentalia`, 
  `Boreoeutheria`, `Galloanserae`).
  
# taxodist 0.3.0

## New features
* `taxo_search()`: Added a new interactive search function that queries The Taxonomicon and         returns a tidy data frame of all available IDs, ranks, and authors for a given taxon name.
* Direct ID Support: `get_lineage()` and all distance functions (e.g.,                             `taxo_distance()`, `distance_matrix()`) now accept Taxonomicon numeric IDs. This                 provides a fallback for computing distances when homonyms or historical ranks cannot be          disambiguated by name alone.

## Network Handling
* Added a 30-second timeout to all HTTP requests to prevent the R session from hanging             indefinitely when The Taxonomicon servers are overloaded.
* Network failures, offline servers, or bad HTTP statuses now immediately emit a clear,            informative warning to the user rather than failing silently.

## Bug fixes
* Extended rank-prefix filter to remove additional bare rank tokens that were
  being retained as spurious lineage nodes: `Subgenus`, `Section`, `Division`,
  `Subdivision`, `Supercohort`, `Infracohort`, `Subsection`, `Candidatus`,
  `Parvphylum`, `Branch`, and `Go to` (a navigation artefact from The
  Taxonomicon page layout).
* Duplicate lineage nodes originating from data-quality issues in The
  Taxonomicon (e.g. `Uropygi` appearing twice in Thelyphonida, `Myxomycetes`
  appearing twice in Physarum) are now collapsed to a single occurrence.
  A warning is emitted when deduplication occurs so the user is aware of the
  upstream data issue.
* `get_taxonomicon_id()` now collects all biological matches for a taxon name
  before returning, and emits a warning when multiple valid biological entries
  are found in The Taxonomicon. This surfaces homonym ambiguity (e.g. `Nereis`
  matching both a polychaete worm and a butterfly genus) that was previously
  silent.
* `get_taxonomicon_id()` now properly follows taxonomic redirects and removes interface            noise (`N|T|P...`) from warnings.
* When multiple biological entries are found (e.g. `Bacteria`), the warning now lists all          available numeric IDs without duplicates, allowing users to make an informed choice.
* Made lineage parsing more robust by strictly truncating any philosophical or pre-basal           nodes (e.g. `organisms`) appearing before `Biota`.
* Added validation to `get_lineage_by_id()` to silently return `NULL` for non-numeric              strings, preventing upstream server fallback errors.
* Fixed a lineage parsing issue where the parent `Genus` was erroneously truncated from            species-level queries. Distances and Most Recent Common Ancestors (MRCAs) between congeneric     species are now computed accurately.
* Fixed an issue where taxonomic author names enclosed in double quotes (e.g., `"Redtenbacher,     1906"`) were evading the lineage cleaning pipeline and being erroneously retained as part of     the clade name.

# taxodist 0.4.0

## New features

* `cache_info()`: inspects the current session cache, reporting the number of
  cached lineages and IDs, total memory used, and the names of all taxa whose
  lineages are stored. Returns a named list invisibly for programmatic use.
* `focal_distances()`: computes taxonomic distances from a single focal taxon
  to all members of a community vector. Returns a sorted data frame with
  columns `taxon`, `distance`, `mrca`, and `mrca_depth`, making it easy to
  identify the closest and most distant relatives of a focal species in a
  community context.

## Bug fixes

* Fixed a lineage parsing bug where taxa with `Subphylum` and `Infraphylum`
  rank prefixes (e.g., Craniata, Vertebrata) were being incorrectly discarded.
  The dagger symbol (`†`) is now stripped before rank prefixes are removed,
  ensuring that prefixed names like `Subphylum † Craniata` are correctly
  parsed to `Craniata` rather than falling through as bare rank tokens.
  `Subphylum` and `Infraphylum` have been removed from `bare_ranks` accordingly.
* Fixed parsing of `auct.` author annotations (e.g., `Subphylum Craniata auct.`)
  that were preventing rank-prefixed names from being correctly cleaned.
* Fixed intermediate clades being silently dropped from lineages when their
  author strings contained ampersands, hyphens (e.g., Cavalier-Smith), or
  trailing commas without a year (e.g., `Romeriida Gauthier,`). The author
  removal pipeline now handles these patterns correctly.
* Fixed a bug where the dagger extinction marker (`†`) appearing between a
  rank prefix and a taxon name (e.g., `Family † Dromaeosauridae`) caused the
  taxon to be lost entirely. Dagger removal now precedes all other cleaning
  steps.
  
# taxodist 0.5.0

## New features

* `taxobase`: Added a built-in reference dataset containing 50 taxonomic clades. 
   The dataset includes pre-computed distance matrices, lineage paths, search 
   outputs, and coverage vectors. It acts as a complete offline fallback, allowing 
   package examples, vignettes, and unit tests to run instantly without requiring 
   a network connection to The Taxonomicon.

## Bug fixes

  * Fixed a bug where any taxon whose name contains the substring "planet" (e.g.               *Periplaneta*) was incorrectly discarded by the astronomical homonym filter in               `get_taxonomicon_id()` and `taxo_search()`. The filter now uses word boundaries              (`\bplanet\b`) so that genus names containing "planet" as an infix are no longer mistakenly   treated as non-biological entries.
  * Fixed a lineage parsing bug affecting genera whose Taxonomicon page does not link back to   the genus itself (e.g. *Periplaneta*, which only links to its constituent species). The      lineage is now correctly truncated before child species entries, and the genus name is       appended if absent.
  * Added `"Epifamily"` to the bare rank token filter so that unprefixed rank labels are not   retained as spurious lineage nodes.
  * Fixed a lineage parsing bug where the superscript type-marker `ᵀ` (U+1D40) immediately     following a taxon name prevented word-boundary matching in the cutoff search, causing the    lineage to overshoot the target node and incorrectly include descendant taxa as ancestors.   This silently affected any type genus or type species (i.e. the majority of genera), and     was particularly visible in deep clades such as *Dinosauria* where several child clades      were being parsed as ancestor nodes. Higher-rank taxa (clades, orders, families — anything   above genus level) were disproportionately affected.
  * Fixed a lineage parsing bug where navigation header text prepended to the                  `#divPageContent` block could contain the target taxon name (e.g. in a cited reference       title), causing the lineage cutoff to land in the junk header rather than the actual         taxonomic tree. Affected taxa received a garbage two- or three-line lineage and were         silently discarded by the `"Biota"` membership check in `get_taxonomicon_id()`.              *Drosophila* (the fly genus, ID 28940) was one such case.

# taxodist 0.6.0

## Breaking changes

* The distance definition is now a proper ultrametric. Distance is zero only
  when two complete lineages identify the same hierarchy node. Distinct
  ancestor-descendant pairs now receive `1 / depth(MRCA)` instead of zero.
  Use `is_member()` or `taxo_path()` when the question concerns clade
  membership or containment rather than distance between nodes.
* Package terminology now consistently describes the measure as a taxonomic
  hierarchy distance rather than a phylogenetic distance.

## Data and documentation

* `taxobase` is now distributed as a single documented list with provenance
  metadata, a compact reference distance matrix, and an offline matrix used to
  build the statistical applications vignette reproducibly.
* Added a methodological vignette defining the distance, demonstrating its
  ultrametric property, and documenting its assumptions and limitations.
* Revised the introduction and statistical applications vignettes to remove
  stale hand-written numerical output and distinguish taxonomic dendrograms
  from inferred phylogenetic trees.
* Added a package citation for taxodist and a separate citation for The
  Taxonomicon data source.

## Bug fixes and robustness

* MRCA detection now uses only the continuous common lineage prefix, preventing
  repeated names after divergence from being mistaken for shared ancestry.
* `is_member()` now performs exact, case-insensitive matching and treats regular
  expression characters literally.
* `taxo_path()` now correctly handles identical taxa, ancestor-descendant pairs,
  path direction, lineage depths, and pairs with no common ancestor.
* Empty and single-taxon inputs are handled safely by distance matrices,
  clustering, ordination, heatmaps, closest-relative queries, and focal
  distances. Ordination also reduces `k` when required by sample size.
* Clustering, ordination, and heatmaps now detect `NA` and infinite distances
  and return informative fallback objects instead of failing downstream.
* `load_cache()` validates cache structure, names, and value types before
  modifying the active cache.

# taxodist 0.7.0

## Performance

* Optimized `distance_matrix()` to fill the compact lower-triangular `dist`
  representation directly instead of allocating a full square matrix.
* Added a numeric-only pairwise distance path for matrix construction,
  avoiding allocation of a detailed S3 result object for every taxon pair.
  Public return values and distance semantics are unchanged.
* Cached the final name-resolved lineage used by `distance_matrix()`, avoiding
  repeated filtering and regular-expression work in warm matrix construction.

## Compatibility and validation

* Preserved the public API and the continuous-prefix MRCA distance semantics.
* Revalidated the reference matrices and lineage outputs against the Python
  and Julia implementations.
* Expanded regression coverage for the optimized numeric and cache paths; the
  package test suite now reports complete source coverage.
* Cache inspection and clearing now address the private cache environment
  explicitly, avoiding ambiguous search-list lookup during installed use.

# taxodist 0.8.0

## Auditable name resolution

* Added `taxo_resolve()` for batch name resolution with explicit `resolved`,
  `ambiguous`, `unresolved`, and `retrieval_error` statuses. A missing name is
  therefore distinguishable from a failed network request or lineage fetch.
  The result preserves all candidates considered, the selected Taxonomicon ID,
  the retrieved lineage, and lineage depth for every input.
* Ambiguity handling is now an explicit policy: warn and select the first
  candidate, select it silently, or stop and require numeric IDs.
* `distance_matrix()` now accepts a `taxodist_resolution` object directly,
  reuses its stored lineages without additional network requests, and keeps
  the original queried names as matrix labels.
* Added `taxo_from_lineages()` for curated, unpublished, frozen, or fully
  offline classifications. User-supplied lineages receive deterministic identifiers
  and flow through the same matrix and bundle APIs as online resolutions.

## Portable analysis bundles

* Added `taxo_bundle()` to combine resolution records, retrieved lineages, the
  distance matrix, metric definition, source details, and software metadata in
  one self-contained object.
* Added `write_taxodist_bundle()` and `read_taxodist_bundle()` for a documented
  JSON format designed for exchange among the R, Python, and Julia
  implementations.
* Added a bundled JSON Schema for format version 1.0. Missing distances are
  represented as JSON `null`, while disconnected infinite distances use the
  explicit strings `"Infinity"` and `"-Infinity"`.
* `distance_matrix()` now also accepts a `taxodist_bundle` and returns its
  stored matrix without network access or recomputation.

## Robustness and documentation

* Repeated names in `taxo_resolve()` are now retrieved only once while their
  original order and repetition are preserved in the audit table.
* Bundle validation now checks resolution fields, candidate counts, selected
  IDs, lineage depths, matrix labels, and agreement between stored distances
  and the stored lineages.
* Bundle provenance now supports sources without a URL and preserves custom
  lineage sources through JSON round trips.
* HTTP requests now report the installed package version in their user agent,
  and malformed lineage pages fail cleanly instead of escaping as parser
  errors.
* Reference-data tests now preserve the package version that originally built
  `taxobase` instead of requiring it to equal the current development version.
* Corrected stale generated documentation that still used phylogenetic
  terminology or the pre-0.6 ancestor--descendant distance behavior.

# taxodist 0.9.0

## Bug fixes

* `taxo_resolve()` now truncates each selected lineage at the queried taxon
  exactly as `get_lineage()` does. Previously it stored the raw lineage of the
  Taxonomicon page, so for genera whose page does not list the genus itself
  (e.g. *Periplaneta*) the lineage stopped at the family, `resolved_name`
  reported the family, and `distance_matrix(taxo_resolve(x))` could differ
  from `distance_matrix(x)`. Bundles built from such resolutions inherited the
  difference.
* Cleaned and uncleaned lineages are now cached separately. Previously, a call
  to `get_lineage_by_id(id, clean = FALSE)` was cached under the same key as
  the cleaned lineage, so later calls with `clean = TRUE` returned nodes above
  `Biota` for the rest of the session. A single request now fills both entries.
* Taxon names are normalized before lookup: surrounding whitespace is removed,
  repeated spaces are collapsed, and the queried name is matched
  case-insensitively when truncating the lineage. Previously
  `get_lineage("tyrannosaurus")` or `get_lineage("Tyrannosaurus ")` appended a
  spurious extra node, inflating the lineage depth, and the distance between
  `"Tyrannosaurus"` and `"tyrannosaurus"` was positive.
* `taxo_heatmap()` no longer forwards `...` to both `distance_matrix()` and
  `stats::heatmap()`. Graphical arguments such as `col` previously caused an
  "unused argument" error when taxon names were supplied, and `progress`
  produced a graphical-parameter warning. `verbose` and `progress` are now
  explicit arguments.
* `cache_info()` now reports the names of cached taxa, together with their
  Taxonomicon IDs in a new `ids` element. It previously listed only the IDs
  under `taxa`.
* `get_lineage_by_id()` now gives an informative error for a vector of several
  IDs instead of a low-level coercion error, accepts IDs supplied as numbers,
  and ignores surrounding whitespace. `get_lineage()` likewise requires a
  single name.
* Text scraped from The Taxonomicon is escaped before being included in
  warnings, so entries containing braces can no longer break the homonym
  warning.
* Search URLs now percent-encode reserved characters such as `&`.
* `focal_distances()` reports the terminal lineage node as the MRCA of a
  community member identical to the focal taxon, consistently with
  `taxo_distance()`, rather than echoing the input name.
* `closest_relative()` returns consecutive row names after sorting.

## Behavior changes

* `check_coverage()` returns `NA` when the answer is unknown because The
  Taxonomicon could not be reached or a candidate lineage could not be
  retrieved. Previously such names were reported as `FALSE`, i.e. as absent
  from the database.
* Validating a bundle when it is printed or passed to `distance_matrix()` no
  longer recomputes every distance from the stored lineages; this full
  integrity check is performed by `write_taxodist_bundle()` and
  `read_taxodist_bundle()`. Printing a bundle of 800 taxa went from about 5 s
  to under 0.2 s.

## Performance

* `distance_matrix()` now computes distances in vectorized form. Every lineage
  node is mapped to an integer identifying its complete path from the root, so
  the MRCA depth of two taxa is the number of depths at which these
  identifiers agree. For 800 taxa the computation dropped from about 14 s
  (with the default progress bar) to 0.3 s, with identical results.
* The progress bar of `distance_matrix()` is updated once per column instead
  of once per pair; previously it tripled the running time.
* HTTP requests are retried with exponential backoff after network errors and
  server errors (5xx). The number of attempts is controlled by the new
  `taxodist.max_tries` option (default 3).

## New features

* `taxo_tree()` builds an ape-compatible `phylo` object directly from lineages.
  Its internal nodes are the named clades at which the taxa diverge, and its
  branch lengths make the tree ultrametric with patristic distances equal to
  `distance_matrix()`.
* `cross_distances()` computes the rectangular matrix of distances between
  two sets of taxa without computing the full square matrix.
* `compare_bundles()` reports which taxa changed status, identifier, or
  lineage between two bundles, and which pairwise distances changed, e.g. when
  an analysis is repeated after The Taxonomicon is updated. Bundle files can
  be passed directly by path.
* `taxo_resolve()` and `taxo_bundle()` gain a `within` argument that restricts
  candidates to a clade (for example `within = "Animalia"`), resolving many
  cross-kingdom homonyms automatically.
* `save_cache()` and `load_cache()` default to a persistent per-user file in
  `tools::R_user_dir("taxodist", "cache")`. `save_cache()` now returns the
  written path invisibly.
* `distance_matrix()`, `cross_distances()`, and `taxo_tree()` accept a
  `taxodist_bundle` or `taxodist_resolution` wherever lineages are needed.

## Internal changes

* Removed the dependency on purrr and replaced superseded rvest functions
  (`html_nodes()`, `html_node()`) with `html_elements()` and `html_element()`.
  The unused suggested package webmockr was also removed.
* The test suite is split into one file per area (retrieval, resolution,
  distances, trees, bundles, cache, analysis, and live network tests).
* The parsing of Taxonomicon search pages is shared by
  `get_taxonomicon_id()`, `taxo_search()`, and `taxo_resolve()`.
* The package citation now takes its version from the installed package
  metadata.
* Tests were extended to cover every new behavior and the package retains
  complete source coverage.

# taxodist 1.0.0

First stable release. The default distance metric, the other metrics, and the
bundle format (schema 1.1) are now stable: later 1.x releases will keep their
definitions and continue to read bundles written by this version.

## New features

* **Ranks.** The rank of every lineage node (e.g. `"Kingdom"`, `"Clade"`,
  `"Superfamily"`, `"Genus"`) is now read from The Taxonomicon and cached
  with the lineage. The new `lineage_ranks()` returns each node with its rank
  and depth. `taxo_resolve()` results gain a `ranks` list column, and
  `taxo_from_lineages()` accepts user-supplied ranks through a new `ranks`
  argument.
* **Distance metrics.** `taxo_distance()`, `distance_matrix()`,
  `cross_distances()`, `closest_relative()`, `focal_distances()`,
  `taxo_heatmap()`, and `taxo_bundle()` gain a `metric` argument, also
  accepted by `taxo_cluster()` and `taxo_ordinate()` through `...`:
  * `"reciprocal"` (default): `1 / depth(MRCA)`, unchanged;
  * `"path"`: the number of edges between the two nodes in the
    classification tree;
  * `"wu_palmer"`: `1 - 2 depth(MRCA) / (depth(A) + depth(B))`, which
    normalizes for the different resolution of different parts of the
    classification;
  * `"ranked"`: half the path length through the standard Linnaean ranks
    only (kingdom to species, configurable with the `taxodist.rank_levels`
    option), i.e. the step weights of Clarke and Warwick.

  Distance matrices record the metric in their `"method"` attribute, which
  `stats::hclust()` reports.
* **Taxonomic distinctness.** `taxo_distinctness()` computes Warwick and
  Clarke's taxonomic diversity (Δ), distinctness (Δ*), average taxonomic
  distinctness (Δ+), and its variation (Λ+) for one or more samples, from
  presence or abundance data. Results agree with `vegan::taxondive()`.
  `taxo_distinctness_test()` compares Δ+ and Λ+ of a sample with random
  samples of the same size from a master list.
* `taxonomicon_available()` reports whether The Taxonomicon can be reached,
  with a short timeout. Examples that need the server now run only when it
  answers (`@examplesIf`), so checking the package no longer waits on network
  timeouts when the server is down.
* **Bundle schema 1.1.** Bundles record the metric used (and, for
  `"ranked"`, its rank levels) and the lineage ranks. Distances are verified
  with the stored metric. Schema 1.0 bundles are still read.

## Breaking changes: lineage parsing

* Lineages are now read from the structure of The Taxonomicon page instead of
  its text. Each node's name is taken from the page's own markup for the
  scientific name, so author strings, rank labels, extinction and type
  markers, and vernacular names can no longer leak into, or remove, nodes.
  Previously about 16% of lineages contained malformed nodes such as
  `"Passeres (C."`, `"Supramyomorpha D'Elía,"`, `"Delphinida de"`,
  `"Alismatales R.Br. ex"`, or `"Bryophyta A.Braun"`; because these contained
  spaces, name-based queries silently dropped them, so lineage depths, and
  therefore distances, differed from the classification and from ID-based
  queries. Some clades were also lost entirely (for example
  `Terrestrornithes`, `Dendrornithes`, and `Anomalogonates` in passerine
  birds), taxa whose page row is labelled `[plesion]` or `Anamorph genus`
  were recorded under that label, and the `‡` marker was not removed.
  **Distances computed with this version can differ from earlier versions.**
* Species are now recorded as nodes of their subspecies, and subgenera as
  nodes of their species, when they appear in the lineage. Previously
  `"Subspecies Homo sapiens sapiens"` was parsed as `"Subspecies"` and the
  species node was dropped.
* Pages without the structured hierarchy block are still parsed from text,
  now by extracting the scientific name at the start of each line rather than
  by deleting author strings, which handles lower-case particles (`de`,
  `van der`, `d'`), non-ASCII initials (`Özdikmen`), and abbreviations
  (`R.Br. ex`).
* The `taxobase` reference dataset was rebuilt with the new parser.

## Bug fixes

* Taxon names containing regular-expression metacharacters, notably
  parentheses in subgenus names such as `"Microtus (Agricola)"`, are now
  matched literally. Whole-word matching no longer fails next to
  punctuation and recognises non-ASCII letters.
* `taxo_resolve()` reports a name whose candidates were all retrieved but are
  not biological (for example astronomical objects) as `"unresolved"`
  instead of `"retrieval_error"`.
* When a name has several biological entries, entries marked as invalid,
  rejected, nude, or illegitimate names (`nom. inval.`, `nom. rej.`,
  `nom. nud.`, `nom. illeg.`) are now considered only after the others.
  Previously `"Amanita"` resolved to an invalid name instead of the
  conserved genus.
* `load_cache()` discards lineages from cache files written by earlier
  versions, whose parser could produce malformed nodes, and keeps their
  taxon IDs. `save_cache()` records the parser version in the file.

## Network

* Consecutive requests to The Taxonomicon are spaced by at least
  `getOption("taxodist.min_interval")` seconds (default 0.25), so large
  batches do not overload the server. Set it to 0 to disable.

## Internal changes

* xml2 moved from Suggests to Imports (it was already required by rvest).
* Trimmed copies of real Taxonomicon pages are included as test fixtures,
  together with their expected lineages. The same fixtures are used by the
  Python and Julia implementations.
