# Session cache: save, load, inspect

test_that("clear_cache returns invisible NULL", {
  expect_invisible(clear_cache())
})

test_that("save_cache creates a file with cache contents", {
  clear_cache()
  assign("id_Carnotaurus", "12345", envir = taxodist:::.taxodist_cache)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  expect_invisible(save_cache(tmp))
  expect_true(file.exists(tmp))
})

test_that("load_cache restores entries into the cache", {
  clear_cache()
  assign("id_Carnotaurus", "12345", envir = taxodist:::.taxodist_cache)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  save_cache(tmp)
  clear_cache()
  expect_false(exists("id_Carnotaurus", envir = taxodist:::.taxodist_cache))
  load_cache(tmp)
  expect_true(exists("id_Carnotaurus", envir = taxodist:::.taxodist_cache))
  expect_equal(get("id_Carnotaurus", envir = taxodist:::.taxodist_cache), "12345")
})

test_that("save_cache / load_cache round-trip preserves all entries", {
  clear_cache()
  assign("id_Tyrannosaurus", "50841", envir = taxodist:::.taxodist_cache)
  assign("lin_50841", c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  save_cache(tmp)
  clear_cache()
  load_cache(tmp)
  expect_equal(get("id_Tyrannosaurus",  envir = taxodist:::.taxodist_cache), "50841")
  expect_equal(get("lin_50841", envir = taxodist:::.taxodist_cache),
               c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus"))
})

test_that("load_cache errors on missing file", {
  expect_error(load_cache("nonexistent_file.rds"))
})

test_that("load_cache rejects a non-list RDS without modifying cache", {
  clear_cache()
  assign("id_existing", "123", envir = taxodist:::.taxodist_cache)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(42, tmp)
  expect_error(load_cache(tmp), "expected a named list")
  expect_equal(
    get("id_existing", envir = taxodist:::.taxodist_cache),
    "123"
  )
})

test_that("load_cache rejects unnamed cache entries", {
  clear_cache()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(list("12345"), tmp)
  expect_error(load_cache(tmp), "unique, non-empty names")
})

test_that("load_cache rejects invalid ID and lineage values", {
  clear_cache()
  tmp_id <- tempfile(fileext = ".rds")
  tmp_lin <- tempfile(fileext = ".rds")
  tmp_matrix_lin <- tempfile(fileext = ".rds")
  on.exit(unlink(c(tmp_id, tmp_lin, tmp_matrix_lin)))

  saveRDS(list(id_bad = 12345), tmp_id)
  expect_error(load_cache(tmp_id), "ID entries")

  saveRDS(list(lin_bad = 12345), tmp_lin)
  expect_error(load_cache(tmp_lin), "lineage entries")

  saveRDS(list(matrix_lineage_bad = 12345), tmp_matrix_lin)
  expect_error(load_cache(tmp_matrix_lin), "matrix lineage entries")
})

test_that("save_cache returns invisibly", {
  clear_cache()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  expect_invisible(save_cache(tmp))
})

test_that("load_cache returns invisible NULL", {
  clear_cache()
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  save_cache(tmp)
  clear_cache()
  expect_invisible(load_cache(tmp))
  expect_equal(ls(taxodist:::.taxodist_cache), character(0))
})

test_that("cache is used on second call to get_taxonomicon_id", {
  lookup_id <- taxodist:::.lookup_taxonomicon_id
  get_id <- function(...) lookup_id(...)$id
  clear_cache()
  assign("id_Tyrannosaurus", "50841", envir = taxodist:::.taxodist_cache)
  call_count <- 0L
  mockery::stub(lookup_id, "httr::RETRY", function(...) {
    call_count <<- call_count + 1L
    NULL
  })
  result <- get_id("Tyrannosaurus")
  expect_equal(call_count, 0L)
  expect_equal(result, "50841")
})

test_that("get_taxonomicon_id verbose prints messages on cache hit", {
  clear_cache()
  assign("id_Carnotaurus", "99999", envir = taxodist:::.taxodist_cache)
  expect_no_error(get_taxonomicon_id("Carnotaurus", verbose = TRUE))
})

test_that("get_lineage_by_id verbose prints messages on cache hit", {
  clear_cache()
  assign("lin_99999",
         c("Biota", "Animalia", "Dinosauria", "Abelisauridae", "Carnotaurus"),
         envir = taxodist:::.taxodist_cache)
  expect_no_error(get_lineage_by_id("99999", verbose = TRUE))
})

test_that("get_lineage_by_id cache hit with verbose prints message", {
  clear_cache()
  assign("lin_50841",
         c("Biota", "Animalia", "Dinosauria", "Theropoda", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  result <- get_lineage_by_id("50841", verbose = TRUE)
  expect_type(result, "character")
})

test_that("cache_info returns invisible list with correct structure", {
  clear_cache()
  result <- cache_info()
  expect_invisible(cache_info())
  expect_type(result, "list")
  expect_named(result, c("n_lineages", "n_ids", "taxa", "ids", "size_bytes"))
})

test_that("cache_info reports zero counts on empty cache", {
  clear_cache()
  result <- cache_info()
  expect_equal(result$n_lineages, 0L)
  expect_equal(result$n_ids, 0L)
  expect_equal(result$taxa, character(0))
})

test_that("cache_info counts lineages and IDs correctly", {
  clear_cache()
  assign("id_Tyrannosaurus", "50841", envir = taxodist:::.taxodist_cache)
  assign("id_Velociraptor",  "12345", envir = taxodist:::.taxodist_cache)
  assign("lin_50841", c("Biota", "Animalia", "Dinosauria", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  result <- cache_info()
  expect_equal(result$n_lineages, 1L)
  expect_equal(result$n_ids, 2L)
})

test_that("cache_info reports taxon names alongside their IDs", {
  clear_cache()
  assign("lin_50841", c("Biota", "Animalia", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  assign("lin_12345", c("Biota", "Animalia", "Velociraptor"),
         envir = taxodist:::.taxodist_cache)
  assign("lin_999", character(0), envir = taxodist:::.taxodist_cache)
  result <- cache_info()
  expect_equal(result$ids, c("12345", "50841", "999"))
  expect_equal(result$taxa, c("Velociraptor", "Tyrannosaurus", NA_character_))
})

test_that("cache_info size_bytes is numeric and positive when cache is populated", {
  clear_cache()
  assign("lin_50841", c("Biota", "Animalia", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  result <- cache_info()
  expect_type(result$size_bytes, "double")
  expect_gt(result$size_bytes, 0)
})

test_that("cache_info prints without error", {
  clear_cache()
  assign("lin_50841", c("Biota", "Animalia", "Tyrannosaurus"),
         envir = taxodist:::.taxodist_cache)
  expect_no_error(cache_info())
})


test_that("save_cache and load_cache use a persistent per-user file by default", {
  cache_root <- tempfile("taxodist-cache-")
  old <- Sys.getenv("R_USER_CACHE_DIR", unset = NA)
  Sys.setenv(R_USER_CACHE_DIR = cache_root)
  on.exit({
    if (is.na(old)) Sys.unsetenv("R_USER_CACHE_DIR") else Sys.setenv(R_USER_CACHE_DIR = old)
    unlink(cache_root, recursive = TRUE)
    clear_cache()
  }, add = TRUE)

  clear_cache()
  assign("lin_1", c("Biota", "Alpha"), envir = taxodist:::.taxodist_cache)
  path <- save_cache()
  expect_true(file.exists(path))
  expect_true(startsWith(normalizePath(path), normalizePath(cache_root)))

  clear_cache()
  load_cache()
  expect_equal(get("lin_1", envir = taxodist:::.taxodist_cache), c("Biota", "Alpha"))
})

test_that("load_cache validates uncleaned lineage entries", {
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(list(rawlin_1 = 1), tmp)
  expect_error(load_cache(tmp), "lineage entries")
})
