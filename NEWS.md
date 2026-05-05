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
