#' Taxonomic diversity and distinctness indices
#'
#' Computes the taxonomic diversity and distinctness indices of Warwick and
#' Clarke for one or more samples: taxonomic diversity \eqn{\Delta},
#' taxonomic distinctness \eqn{\Delta^*}, average taxonomic distinctness
#' \eqn{\Delta^+}, and variation in taxonomic distinctness \eqn{\Lambda^+}.
#'
#' @param comm The samples. One of:
#'   * a character vector of taxon names (one sample, presence only);
#'   * a named numeric vector of abundances, named by taxon (one sample);
#'   * a matrix or data frame of abundances or presences with one row per
#'     sample and one column per taxon (column names are taxon names).
#' @param d Optional `dist` object whose labels include every taxon present
#'   in `comm`, for example from [distance_matrix()]. If `NULL` (default),
#'   distances are computed with [distance_matrix()] and `metric`.
#' @param metric Distance metric used when `d` is `NULL`; see
#'   [taxo_distance()]. The default `"ranked"` gives the step weights of
#'   Clarke and Warwick (1 for congeneric taxa, 2 for confamilial, ...).
#' @param verbose,progress Passed to [distance_matrix()] when `d` is `NULL`.
#'
#' @return A data frame of class `"taxodist_distinctness"` with one row per
#'   sample and columns:
#' \describe{
#'   \item{`sample`}{Sample name (row name of `comm`, or `"sample"`).}
#'   \item{`S`}{Number of taxa present.}
#'   \item{`Delta`}{Taxonomic diversity \eqn{\Delta}.}
#'   \item{`Delta_star`}{Taxonomic distinctness \eqn{\Delta^*}.}
#'   \item{`Delta_plus`}{Average taxonomic distinctness \eqn{\Delta^+}.}
#'   \item{`Lambda_plus`}{Variation in taxonomic distinctness
#'     \eqn{\Lambda^+}.}
#'   \item{`SDelta_plus`}{Total taxonomic distinctness \eqn{S\Delta^+}.}
#' }
#' The indices are `NA` for samples with fewer than two taxa, or when a
#' distance between present taxa is missing or infinite.
#'
#' @details
#' With \eqn{x_i} the abundance of taxon \eqn{i}, \eqn{n = \sum x_i},
#' \eqn{S} the number of taxa present, and \eqn{\omega_{ij}} the distance
#' between taxa \eqn{i} and \eqn{j}, sums running over pairs \eqn{i < j}:
#' \deqn{\Delta = \sum\sum \omega_{ij} x_i x_j / [n(n-1)/2]}
#' \deqn{\Delta^* = \sum\sum \omega_{ij} x_i x_j / \sum\sum x_i x_j}
#' \deqn{\Delta^+ = \sum\sum \omega_{ij} / [S(S-1)/2]}
#' \deqn{\Lambda^+ = \sum\sum (\omega_{ij} - \Delta^+)^2 / [S(S-1)/2]}
#' \eqn{\Delta^+} and \eqn{\Lambda^+} use presence only and, unlike
#' species richness, their expected values do not depend on sampling
#' effort; see [taxo_distinctness_test()]. With presence-only data,
#' \eqn{\Delta} and \eqn{\Delta^*} are reported for completeness and equal
#' \eqn{\Delta^+}.
#'
#' Any metric can be used. The `"ranked"` metric reproduces the classical
#' indices, which use equal steps between Linnaean ranks; the other metrics
#' use the full resolution of The Taxonomicon classification.
#'
#' @references
#' Warwick, R.M. & Clarke, K.R. (1995). New 'biodiversity' measures reveal a
#' decrease in taxonomic distinctness with increasing stress. *Marine
#' Ecology Progress Series* 129, 301--305.
#'
#' Clarke, K.R. & Warwick, R.M. (1998). A taxonomic distinctness index and
#' its statistical properties. *Journal of Applied Ecology* 35, 523--531.
#'
#' Clarke, K.R. & Warwick, R.M. (2001). A further biodiversity index
#' applicable to species lists: variation in taxonomic distinctness. *Marine
#' Ecology Progress Series* 216, 265--278.
#'
#' @seealso [taxo_distinctness_test()], [distance_matrix()]
#' @export
#' @examples
#' lineages <- list(
#'   A = c("Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "A"),
#'   B = c("Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "B"),
#'   C = c("Biota", "Animalia", "Chordata", "Mammalia", "Primates", "C"),
#'   D = c("Biota", "Animalia", "Arthropoda", "Insecta", "Diptera", "D")
#' )
#' d <- distance_matrix(taxo_from_lineages(lineages), metric = "path")
#' comm <- rbind(
#'   site1 = c(A = 5, B = 3, C = 1, D = 0),
#'   site2 = c(A = 1, B = 0, C = 2, D = 4)
#' )
#' taxo_distinctness(comm, d)
taxo_distinctness <- function(comm, d = NULL, metric = "ranked",
                              verbose = FALSE, progress = TRUE) {
  samples <- .distinctness_samples(comm)
  taxa <- colnames(samples)
  if (is.null(d)) {
    metric <- .match_metric(metric)
    d <- distance_matrix(taxa, verbose = verbose, progress = progress, metric = metric)
  }
  full <- .distinctness_matrix(d, taxa)
  metric_used <- attr(d, "method", exact = TRUE)

  rows <- lapply(seq_len(nrow(samples)), function(k) {
    .distinctness_indices(samples[k, ], full)
  })
  out <- data.frame(
    sample = rownames(samples),
    do.call(rbind, rows),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  structure(
    out,
    class = c("taxodist_distinctness", "data.frame"),
    metric = if (is.null(metric_used)) NA_character_ else metric_used
  )
}

#' Test taxonomic distinctness against random samples from a species pool
#'
#' Compares the average taxonomic distinctness \eqn{\Delta^+} and its
#' variation \eqn{\Lambda^+} of a sample with their distribution in random
#' samples of the same number of taxa drawn from a master list (species
#' pool), as proposed by Clarke and Warwick.
#'
#' @param sample A character vector of the taxa present in the sample.
#' @param d A `dist` object over the master list, for example from
#'   [distance_matrix()]. Every taxon in `sample` must be one of its labels.
#' @param n_iter Number of random samples. Default `999`.
#'
#' @return An object of class `"taxodist_distinctness_test"`, a list with:
#' \describe{
#'   \item{`statistics`}{A data frame with one row for `Delta_plus` and one
#'     for `Lambda_plus`, and columns `observed`, `expected`, `lower`,
#'     `upper` (the 2.5% and 97.5% quantiles of the random samples),
#'     `p_lower`, and `p_upper` (one-sided permutation p-values for values as
#'     low or as high as observed).}
#'   \item{`simulated`}{A matrix of the simulated values, one row per random
#'     sample.}
#'   \item{`S`, `pool_size`, `n_iter`, `metric`}{Test settings.}
#' }
#'
#' @details
#' The expected \eqn{\Delta^+} is the mean distance over all pairs of the
#' master list, which does not depend on the sample size. The expected
#' \eqn{\Lambda^+} is the mean of the simulated values. A value of
#' \eqn{\Delta^+} below the lower limit indicates that the sample is less
#' taxonomically distinct than expected from the pool, a pattern associated
#' with environmental degradation. Use [set.seed()] for reproducible results.
#'
#' @inherit taxo_distinctness references
#' @seealso [taxo_distinctness()]
#' @export
#' @examples
#' lineages <- list(
#'   A = c("Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "A"),
#'   B = c("Biota", "Animalia", "Chordata", "Mammalia", "Carnivora", "B"),
#'   C = c("Biota", "Animalia", "Chordata", "Mammalia", "Primates", "C"),
#'   D = c("Biota", "Animalia", "Arthropoda", "Insecta", "Diptera", "D"),
#'   E = c("Biota", "Animalia", "Arthropoda", "Insecta", "Coleoptera", "E"),
#'   F = c("Biota", "Plantae", "Tracheophyta", "Magnoliopsida", "Fagales", "F")
#' )
#' pool <- distance_matrix(taxo_from_lineages(lineages), metric = "path")
#' set.seed(1)
#' taxo_distinctness_test(c("A", "B", "C"), pool, n_iter = 199)
taxo_distinctness_test <- function(sample, d, n_iter = 999) {
  if (!inherits(d, "dist")) {
    cli::cli_abort("{.arg d} must be a {.cls dist} object over the master list.")
  }
  if (!is.character(sample) || anyNA(sample)) {
    cli::cli_abort("{.arg sample} must be a character vector of taxon names.")
  }
  if (!is.numeric(n_iter) || length(n_iter) != 1L || is.na(n_iter) ||
      n_iter < 1 || n_iter != floor(n_iter)) {
    cli::cli_abort("{.arg n_iter} must be a single positive integer.")
  }
  sample <- unique(sample)
  full <- as.matrix(d)
  pool <- rownames(full)
  missing <- setdiff(sample, pool)
  if (length(missing) > 0L) {
    cli::cli_abort(c(
      "x" = "Every sampled taxon must belong to the master list.",
      "i" = "Not in {.arg d}: {.val {missing}}."
    ))
  }
  if (any(!is.finite(full[lower.tri(full)]))) {
    cli::cli_abort("The master list distances must all be finite.")
  }
  S <- length(sample)
  if (S < 2L || S > length(pool)) {
    cli::cli_abort(
      "The sample must contain between 2 and {length(pool)} taxa of the master list."
    )
  }

  stats_of <- function(idx) {
    sub <- full[idx, idx, drop = FALSE]
    omega <- sub[lower.tri(sub)]
    delta_plus <- mean(omega)
    c(Delta_plus = delta_plus, Lambda_plus = mean((omega - delta_plus)^2))
  }
  observed <- stats_of(match(sample, pool))
  simulated <- t(vapply(seq_len(n_iter), function(i) {
    stats_of(sample.int(length(pool), S))
  }, numeric(2)))
  colnames(simulated) <- names(observed)

  expected <- c(
    Delta_plus = mean(full[lower.tri(full)]),
    Lambda_plus = mean(simulated[, "Lambda_plus"])
  )
  tol <- sqrt(.Machine$double.eps)
  statistics <- data.frame(
    statistic = names(observed),
    observed = unname(observed),
    expected = unname(expected),
    lower = unname(apply(simulated, 2L, stats::quantile, probs = 0.025)),
    upper = unname(apply(simulated, 2L, stats::quantile, probs = 0.975)),
    p_lower = unname(vapply(names(observed), function(s) {
      (1 + sum(simulated[, s] <= observed[[s]] + tol)) / (n_iter + 1)
    }, numeric(1))),
    p_upper = unname(vapply(names(observed), function(s) {
      (1 + sum(simulated[, s] >= observed[[s]] - tol)) / (n_iter + 1)
    }, numeric(1))),
    stringsAsFactors = FALSE
  )
  metric <- attr(d, "method", exact = TRUE)
  structure(
    list(
      statistics = statistics,
      simulated = simulated,
      S = S,
      pool_size = length(pool),
      n_iter = as.integer(n_iter),
      metric = if (is.null(metric)) NA_character_ else metric
    ),
    class = "taxodist_distinctness_test"
  )
}

#' Print taxonomic distinctness results
#'
#' @param x A `taxodist_distinctness` object from [taxo_distinctness()] or a
#'   `taxodist_distinctness_test` object from [taxo_distinctness_test()].
#' @param ... Additional arguments passed to [print.data.frame()].
#' @return Invisibly returns `x`.
#' @export
print.taxodist_distinctness <- function(x, ...) {
  metric <- attr(x, "metric", exact = TRUE)
  cli::cli_h2("Taxonomic distinctness")
  if (!is.null(metric) && !is.na(metric)) cli::cli_text("Metric: {metric}")
  print.data.frame(x, row.names = FALSE, ...)
  invisible(x)
}

#' @rdname print.taxodist_distinctness
#' @export
print.taxodist_distinctness_test <- function(x, ...) {
  cli::cli_h2("Taxonomic distinctness test")
  cli::cli_text(
    "{x$S} taxa vs {x$n_iter} random samples from a pool of {x$pool_size}",
    if (!is.na(x$metric)) " (metric: {x$metric})" else ""
  )
  print.data.frame(x$statistics, row.names = FALSE, digits = 4, ...)
  invisible(x)
}

# -- Internal helpers ---------------------------------------------------------

# Normalize the supported sample formats to a numeric matrix with one row per
# sample and one column per taxon.
.distinctness_samples <- function(comm) {
  if (is.character(comm)) {
    if (anyNA(comm) || any(!nzchar(trimws(comm)))) {
      cli::cli_abort("{.arg comm} cannot contain missing or empty taxon names.")
    }
    taxa <- unique(comm)
    return(matrix(1, nrow = 1L, ncol = length(taxa),
                  dimnames = list("sample", taxa)))
  }
  if (is.data.frame(comm)) comm <- as.matrix(comm)
  if (is.numeric(comm) && is.null(dim(comm))) {
    if (is.null(names(comm))) {
      cli::cli_abort("A numeric {.arg comm} vector must be named by taxon.")
    }
    comm <- matrix(comm, nrow = 1L, dimnames = list("sample", names(comm)))
  }
  if (!is.matrix(comm) || !(is.numeric(comm) || is.logical(comm))) {
    cli::cli_abort(paste(
      "{.arg comm} must be a character vector of taxa, a named numeric",
      "vector, or a numeric matrix or data frame."
    ))
  }
  storage.mode(comm) <- "double"
  if (is.null(colnames(comm)) || anyNA(colnames(comm)) ||
      anyDuplicated(colnames(comm))) {
    cli::cli_abort("The columns of {.arg comm} must be named by unique taxa.")
  }
  if (anyNA(comm) || any(comm < 0)) {
    cli::cli_abort("Abundances in {.arg comm} must be non-negative and not missing.")
  }
  if (is.null(rownames(comm))) {
    rownames(comm) <- if (nrow(comm) == 1L) "sample" else paste0("sample", seq_len(nrow(comm)))
  }
  comm
}

# Full distance matrix restricted to and ordered by `taxa`.
.distinctness_matrix <- function(d, taxa) {
  if (!inherits(d, "dist")) {
    cli::cli_abort("{.arg d} must be a {.cls dist} object.")
  }
  full <- as.matrix(d)
  missing <- setdiff(taxa, rownames(full))
  if (length(missing) > 0L) {
    cli::cli_abort(c(
      "x" = "{.arg d} has no distances for some taxa.",
      "i" = "Missing: {.val {missing}}."
    ))
  }
  full[taxa, taxa, drop = FALSE]
}

.distinctness_indices <- function(x, full) {
  present <- x > 0
  S <- sum(present)
  empty <- c(S = S, Delta = NA_real_, Delta_star = NA_real_,
             Delta_plus = NA_real_, Lambda_plus = NA_real_,
             SDelta_plus = NA_real_)
  if (S < 2L) return(empty)
  w <- full[present, present, drop = FALSE]
  a <- x[present]
  pairs <- lower.tri(w)
  omega <- w[pairs]
  if (anyNA(omega) || any(!is.finite(omega))) return(empty)

  xx <- outer(a, a)[pairs]
  n <- sum(a)
  delta_plus <- mean(omega)
  c(
    S = S,
    Delta = sum(omega * xx) / (n * (n - 1) / 2),
    Delta_star = sum(omega * xx) / sum(xx),
    Delta_plus = delta_plus,
    Lambda_plus = mean((omega - delta_plus)^2),
    SDelta_plus = S * delta_plus
  )
}
