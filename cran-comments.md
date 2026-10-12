## Submission

This is an update from taxodist 0.8.0 to 0.9.0. It fixes several bugs (most
importantly an inconsistency between name-based and resolution-based lineages),
vectorizes the distance-matrix computation, and adds `taxo_tree()`,
`cross_distances()`, and `compare_bundles()`. See NEWS.md for details.

## R CMD check results
0 errors | 0 warnings | 0 notes

Tested on Windows, Ubuntu, and macOS via GitHub Actions across
R release, devel, oldrel, and oldrel-1.

## Notes
This package requires internet access to retrieve lineage data
from The Taxonomicon (taxonomy.nl). Network tests are skipped
on CRAN and continuous-integration services via `skip_on_cran()` and
`skip_on_ci()`. Unit tests use mocked responses and do not require network
access.

`save_cache()` and `load_cache()` use `tools::R_user_dir("taxodist", "cache")`
only when called without a file argument, i.e. on explicit user request.
Examples, vignettes, and tests never write to the user's home filespace; the
corresponding test redirects `R_USER_CACHE_DIR` to a temporary directory.

The maintainer of The Taxonomicon was contacted to inform them of this 
package and request permission. No response has been received at time 
of submission.
