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
