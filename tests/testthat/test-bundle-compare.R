# ── Bundle validation cost and compare_bundles ──────────────────────────────

make_bundle <- function(lineages, ...) {
  taxo_bundle(taxo_from_lineages(lineages, ...), progress = FALSE)
}

test_that("using a bundle does not recompute its distances", {
  bundle <- make_bundle(list(
    Alpha = c("Biota", "Alpha"),
    Beta = c("Biota", "Beta")
  ))
  # A tampered matrix is still returned as stored by distance_matrix() and
  # print(); the full integrity check runs when writing or reading.
  bundle$matrix[] <- 0.75
  expect_equal(c(distance_matrix(bundle)), 0.75)
  expect_invisible(print(bundle))
  expect_error(
    write_taxodist_bundle(bundle, tempfile(fileext = ".json")),
    "stored distances"
  )
})

test_that("compare_bundles reports unchanged bundles", {
  lineages <- list(Alpha = c("Biota", "A", "Alpha"), Beta = c("Biota", "B", "Beta"))
  comparison <- compare_bundles(make_bundle(lineages), make_bundle(lineages))
  expect_s3_class(comparison, "taxodist_bundle_comparison")
  expect_equal(comparison$taxa$change, c("unchanged", "unchanged"))
  expect_equal(nrow(comparison$distances), 0L)
  expect_invisible(print(comparison))
})

test_that("compare_bundles classifies every kind of change", {
  old <- make_bundle(
    list(
      Same = c("Biota", "Animalia", "Same"),
      Moved = c("Biota", "Animalia", "Moved"),
      Renumbered = c("Biota", "Animalia", "Renumbered"),
      Dropped = c("Biota", "Plantae", "Dropped")
    ),
    ids = c("1", "2", "3", "4")
  )
  new <- make_bundle(
    list(
      Same = c("Biota", "Animalia", "Same"),
      Moved = c("Biota", "Plantae", "Moved"),
      Renumbered = c("Biota", "Animalia", "Renumbered"),
      Added = c("Biota", "Fungi", "Added")
    ),
    ids = c("1", "2", "33", "5")
  )
  comparison <- compare_bundles(old, new)
  changes <- stats::setNames(comparison$taxa$change, comparison$taxa$input)
  expect_equal(
    changes,
    c(Same = "unchanged", Moved = "lineage_changed",
      Renumbered = "id_changed", Dropped = "removed", Added = "added")
  )
  expect_true(is.na(comparison$taxa$status_new[comparison$taxa$input == "Dropped"]))
  expect_equal(comparison$taxa$depth_old[comparison$taxa$input == "Added"], NA_integer_)

  expect_equal(nrow(comparison$distances), 2L)
  expect_setequal(
    paste(comparison$distances$taxon_a, comparison$distances$taxon_b),
    c("Same Moved", "Moved Renumbered")
  )
  expect_true(all(comparison$distances$distance_old == 0.5))
  expect_true(all(comparison$distances$distance_new == 1))
  expect_invisible(print(comparison))
})

test_that("compare_bundles detects status changes and reads bundle files", {
  resolved <- make_bundle(list(Alpha = c("Biota", "Alpha"), Beta = c("Biota", "Beta")))
  failed <- resolved$resolution
  failed$status[[2]] <- "retrieval_error"
  failed$id[[2]] <- NA_character_
  failed$resolved_name[[2]] <- NA_character_
  failed$lineage[2] <- list(NULL)
  failed$lineage_depth[[2]] <- NA_integer_
  failed_bundle <- taxo_bundle(failed, progress = FALSE)

  file <- tempfile(fileext = ".json")
  on.exit(unlink(file), add = TRUE)
  write_taxodist_bundle(failed_bundle, file)

  comparison <- compare_bundles(resolved, file)
  expect_equal(comparison$taxa$change, c("unchanged", "status_changed"))
  expect_equal(comparison$distances$distance_old, 1)
  expect_true(is.na(comparison$distances$distance_new))

  same_missing <- compare_bundles(failed_bundle, failed_bundle)
  expect_equal(nrow(same_missing$distances), 0L)
})

test_that("compare_bundles validates its inputs", {
  bundle <- make_bundle(list(Alpha = c("Biota", "Alpha")))
  expect_error(compare_bundles(list(), bundle), "old")
  expect_error(compare_bundles(bundle, 42), "new")
  single <- compare_bundles(bundle, bundle)
  expect_equal(nrow(single$distances), 0L)
})
