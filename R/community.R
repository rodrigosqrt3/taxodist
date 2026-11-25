#' Taxonomic dissimilarity between samples
#'
#' Computes the taxonomic dissimilarity \eqn{\Gamma^+} or \eqn{\Theta^+} of
#' Clarke, Somerfield and Chapman (2006) between every pair of samples. Unlike
#' Bray--Curtis or Jaccard dissimilarities, they credit "near misses": two
#' samples without any taxon in common are still similar when their taxa are
#' closely related, so assemblages from distant regions, or identified with
#' different taxonomic expertise, can be compared.
#'
#' @inheritParams taxo_distinctness
#' @param comm The samples: a matrix or data frame with one row per sample and
#'   one column per taxon (column names are taxon names). Values greater than
#'   zero mean presence; abundances are not used.
#' @param index `"gamma_plus"` (default) or `"theta_plus"`.
#'
#' @return A `dist` object between samples, with attributes `method` (the
#'   index) and `metric` (the distance metric of `d`). Pairs involving a
#'   sample without taxa, or whose nearest-relative distances are missing or
#'   infinite, are `NA`.
#'
#' @details
#' With \eqn{\omega_{ij}} the distance between taxon \eqn{i} of sample 1 and
#' taxon \eqn{j} of sample 2, and \eqn{s_1}, \eqn{s_2} the numbers of taxa in
#' each sample,
#' \deqn{\Gamma^+ = \left[\sum_i \min_j \omega_{ij} + \sum_j \min_i
#'   \omega_{ij}\right] / (s_1 + s_2)}
#' \deqn{\Theta^+ = \frac{1}{2}\left[\sum_i \min_j \omega_{ij} / s_1 +
#'   \sum_j \min_i \omega_{ij} / s_2\right]}
#' \eqn{\Gamma^+} is the mean distance from every taxon to its nearest
#' relative in the other sample; \eqn{\Theta^+} averages this mean separately
#' within each sample. Taxa shared by both samples contribute zero. When every
#' pair of distinct taxa is at the same distance (a flat classification), the
#' indices reduce to the presence/absence Bray--Curtis (Sorensen) and
#' Kulczynski dissimilarities, times that distance.
#'
#' Clarke et al. scale distances so that the largest possible path length is
#' 100; the indices are then percentages. They are returned here in the units
#' of `d`, so rescale `d` for that convention.
#'
#' @references
#' Clarke, K.R., Somerfield, P.J. & Chapman, M.G. (2006). On resemblance
#' measures for ecological studies, including taxonomic dissimilarities and a
#' zero-adjusted Bray-Curtis coefficient for denuded assemblages. *Journal of
#' Experimental Marine Biology and Ecology* 330, 55--80.
#' \doi{10.1016/j.jembe.2005.12.017}
#'
#' @seealso [taxo_distinctness()], [taxo_ses()]
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
#'   site1 = c(A = 1, B = 0, C = 1, D = 0),
#'   site2 = c(A = 0, B = 1, C = 0, D = 0),
#'   site3 = c(A = 0, B = 0, C = 0, D = 1)
#' )
#' # site1 and site2 share no taxon but are close relatives
#' taxo_dissimilarity(comm, d)
taxo_dissimilarity <- function(comm, d = NULL, index = "gamma_plus",
                               metric = "ranked", verbose = FALSE,
                               progress = TRUE) {
  if (!is.character(index) || length(index) != 1L ||
      !index %in% c("gamma_plus", "theta_plus")) {
    cli::cli_abort("{.arg index} must be {.val gamma_plus} or {.val theta_plus}.")
  }
  samples <- .community_samples(comm)
  full <- .community_distances(samples, d, metric, verbose, progress)
  present <- samples > 0
  n <- nrow(samples)

  pair <- function(k, l) {
    p <- present[k, ]
    q <- present[l, ]
    if (!any(p) || !any(q)) return(NA_real_)
    w <- full[p, q, drop = FALSE]
    to_q <- apply(w, 1L, min)
    to_p <- apply(w, 2L, min)
    if (anyNA(to_q) || anyNA(to_p) || any(!is.finite(c(to_q, to_p)))) {
      return(NA_real_)
    }
    if (identical(index, "gamma_plus")) {
      (sum(to_q) + sum(to_p)) / (length(to_q) + length(to_p))
    } else {
      (mean(to_q) + mean(to_p)) / 2
    }
  }

  values <- numeric(n * (n - 1L) / 2L)
  position <- 0L
  for (l in seq_len(max(n - 1L, 0L))) {
    for (k in (l + 1L):n) {
      position <- position + 1L
      values[[position]] <- pair(k, l)
    }
  }
  structure(
    values,
    Size = n,
    Labels = rownames(samples),
    Diag = FALSE,
    Upper = FALSE,
    method = index,
    metric = attr(full, "metric", exact = TRUE),
    class = "dist"
  )
}

#' Standardized effect sizes of community taxonomic structure
#'
#' Compares the mean pairwise distance (MPD) or mean nearest-taxon distance
#' (MNTD) among the taxa of each sample with its distribution under a null
#' model, giving the standardized effect sizes known as the net relatedness
#' index (NRI) and the nearest taxon index (NTI) of Webb et al. (2002),
#' computed on taxonomic distances.
#'
#' @inheritParams taxo_distinctness
#' @param comm The samples: a character vector of taxa (one sample), a named
#'   numeric vector of abundances, or a matrix or data frame with one row per
#'   sample and one column per taxon.
#' @param d Optional `dist` object whose labels include every taxon of `comm`
#'   and of `pool`. If `NULL`, it is computed with [distance_matrix()].
#' @param index `"mpd"` (default) or `"mntd"`.
#' @param null_model How random communities are generated:
#'   \describe{
#'     \item{`"pool"`}{(default) each sample is replaced by the same number of
#'       taxa drawn at random, without replacement, from `pool`; with
#'       `abundance_weighted = TRUE`, the observed abundances are assigned to
#'       the drawn taxa. This is the null model of Webb et al. (2002).}
#'     \item{`"independent_swap"`}{the presence matrix of all samples is
#'       randomized with `swaps` independent swaps (Gotelli 2000), which keep
#'       the number of taxa in every sample and the number of samples in which
#'       every taxon occurs. Requires at least two samples.}
#'   }
#' @param pool Character vector of taxa forming the regional pool for the
#'   `"pool"` null model. Defaults to every label of `d` (or, when `d` is
#'   `NULL`, every taxon of `comm`).
#' @param n_iter Number of random communities.
#' @param abundance_weighted Logical. If `TRUE`, MPD weights each pair of taxa
#'   by the product of their abundances (including each taxon with itself,
#'   which makes it Rao's quadratic entropy) and MNTD weights each taxon by its
#'   abundance, as in the picante package. Default `FALSE` (presence only).
#' @param swaps Number of successful swaps per random community for
#'   `"independent_swap"`.
#' @param metric Distance metric used when `d` is `NULL`. The default
#'   `"reciprocal"` uses the full depth of the classification.
#'
#' @return A data frame of class `"taxodist_ses"` with one row per sample and
#'   columns `sample`, `S` (number of taxa), `observed`, `null_mean`,
#'   `null_sd`, `z` (`(observed - null_mean) / null_sd`), `NRI` or `NTI`
#'   (`-z`), `p_lower` and `p_upper` (one-sided permutation p-values for
#'   values as low or as high as observed, `(1 + count) / (1 + n_iter)`), and
#'   `n_iter`. Samples with fewer than two taxa have `NA` statistics.
#'
#' @details
#' Positive NRI and NTI (low MPD and MNTD) indicate taxonomic clustering: the
#' taxa of the sample are more closely related than expected from the null
#' model. Negative values indicate overdispersion. NRI reflects the structure
#' of the whole classification and NTI that of its terminal branches. Both
#' depend on the pool, which should be chosen to represent the taxa that could
#' have occurred in the samples. Use [set.seed()] for reproducible results.
#'
#' @references
#' Webb, C.O., Ackerly, D.D., McPeek, M.A. & Donoghue, M.J. (2002).
#' Phylogenies and community ecology. *Annual Review of Ecology and
#' Systematics* 33, 475--505. \doi{10.1146/annurev.ecolsys.33.010802.150448}
#'
#' Gotelli, N.J. (2000). Null model analysis of species co-occurrence
#' patterns. *Ecology* 81, 2606--2621.
#'
#' @seealso [taxo_distinctness_test()], [taxo_dissimilarity()]
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
#' d <- distance_matrix(taxo_from_lineages(lineages), metric = "path")
#' set.seed(1)
#' taxo_ses(c("A", "B", "C"), d, n_iter = 199)
taxo_ses <- function(comm, d = NULL, index = "mpd", null_model = "pool",
                     pool = NULL, n_iter = 999, abundance_weighted = FALSE,
                     swaps = 1000, metric = "reciprocal", verbose = FALSE,
                     progress = TRUE) {
  if (!is.character(index) || length(index) != 1L || !index %in% c("mpd", "mntd")) {
    cli::cli_abort("{.arg index} must be {.val mpd} or {.val mntd}.")
  }
  if (!is.character(null_model) || length(null_model) != 1L ||
      !null_model %in% c("pool", "independent_swap")) {
    cli::cli_abort("{.arg null_model} must be {.val pool} or {.val independent_swap}.")
  }
  .check_count(n_iter, "n_iter")
  .check_count(swaps, "swaps")
  if (!is.logical(abundance_weighted) || length(abundance_weighted) != 1L ||
      is.na(abundance_weighted)) {
    cli::cli_abort("{.arg abundance_weighted} must be TRUE or FALSE.")
  }
  samples <- .community_samples(comm)
  if (is.null(pool)) {
    pool <- if (is.null(d)) colnames(samples) else attr(d, "Labels", exact = TRUE)
  }
  if (!is.character(pool) || anyNA(pool) || length(pool) < 2L) {
    cli::cli_abort("{.arg pool} must be a character vector of at least two taxa.")
  }
  pool <- unique(pool)
  if (length(setdiff(colnames(samples), pool)) > 0L && identical(null_model, "pool")) {
    cli::cli_abort(c(
      "x" = "Every taxon of {.arg comm} must belong to the pool.",
      "i" = "Not in the pool: {.val {setdiff(colnames(samples), pool)}}."
    ))
  }
  taxa <- union(colnames(samples), pool)
  full <- .community_distances(
    matrix(0, 1L, length(taxa), dimnames = list("x", taxa)),
    d, metric, verbose, progress
  )
  if (null_model == "independent_swap" && nrow(samples) < 2L) {
    cli::cli_abort("The {.val independent_swap} null model needs at least two samples.")
  }

  stat <- function(x, labels) {
    .community_statistic(x, full[labels, labels, drop = FALSE], index, abundance_weighted)
  }
  observed <- vapply(seq_len(nrow(samples)), function(k) {
    stat(samples[k, ], colnames(samples))
  }, numeric(1))
  S <- rowSums(samples > 0)

  simulated <- matrix(NA_real_, nrow(samples), n_iter)
  for (it in seq_len(n_iter)) {
    if (null_model == "pool") {
      for (k in seq_len(nrow(samples))) {
        present <- samples[k, ] > 0
        if (sum(present) < 2L) next
        drawn <- pool[sample.int(length(pool), sum(present))]
        simulated[k, it] <- stat(stats::setNames(samples[k, present], drawn), drawn)
      }
    } else {
      random <- .independent_swap(samples, swaps)
      for (k in seq_len(nrow(samples))) {
        simulated[k, it] <- stat(random[k, ], colnames(samples))
      }
    }
  }

  null_mean <- rowMeans(simulated)
  null_sd <- apply(simulated, 1L, stats::sd)
  z <- (observed - null_mean) / null_sd
  z[!is.finite(z)] <- NA_real_
  tol <- sqrt(.Machine$double.eps)
  p_lower <- vapply(seq_along(observed), function(k) {
    if (is.na(observed[[k]])) return(NA_real_)
    (1 + sum(simulated[k, ] <= observed[[k]] + tol)) / (n_iter + 1)
  }, numeric(1))
  p_upper <- vapply(seq_along(observed), function(k) {
    if (is.na(observed[[k]])) return(NA_real_)
    (1 + sum(simulated[k, ] >= observed[[k]] - tol)) / (n_iter + 1)
  }, numeric(1))

  out <- data.frame(
    sample = rownames(samples),
    S = as.integer(S),
    observed = observed,
    null_mean = null_mean,
    null_sd = null_sd,
    z = z,
    stringsAsFactors = FALSE
  )
  out[[if (index == "mpd") "NRI" else "NTI"]] <- -z
  out$p_lower <- p_lower
  out$p_upper <- p_upper
  out$n_iter <- as.integer(n_iter)
  rownames(out) <- NULL
  structure(
    out,
    class = c("taxodist_ses", "data.frame"),
    index = index,
    null_model = null_model,
    abundance_weighted = abundance_weighted,
    pool_size = length(pool),
    metric = attr(full, "metric", exact = TRUE)
  )
}

#' @rdname taxo_ses
#' @param x A `taxodist_ses` object.
#' @param ... Additional arguments passed to [print.data.frame()].
#' @export
print.taxodist_ses <- function(x, ...) {
  index <- attr(x, "index", exact = TRUE)
  cli::cli_h2(
    if (identical(index, "mpd")) "Net relatedness index (MPD)" else "Nearest taxon index (MNTD)"
  )
  metric <- attr(x, "metric", exact = TRUE)
  cli::cli_text(
    "Null model: {attr(x, 'null_model', exact = TRUE)}, pool of ",
    "{attr(x, 'pool_size', exact = TRUE)} taxa",
    if (isTRUE(attr(x, "abundance_weighted", exact = TRUE))) ", abundance-weighted" else "",
    if (!is.null(metric) && !is.na(metric)) "; metric: {metric}" else ""
  )
  print.data.frame(x, row.names = FALSE, digits = 4, ...)
  invisible(x)
}

# -- Internal helpers ---------------------------------------------------------

.check_count <- function(value, arg) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      value < 1 || value != floor(value)) {
    cli::cli_abort("{.arg {arg}} must be a single positive integer.")
  }
}

.community_samples <- function(comm) {
  samples <- .distinctness_samples(comm)
  if (ncol(samples) == 0L) cli::cli_abort("{.arg comm} contains no taxa.")
  samples
}

# Full distance matrix over the taxa of `samples`, computing it when `d` is
# NULL. The metric is kept as an attribute.
.community_distances <- function(samples, d, metric, verbose, progress) {
  taxa <- colnames(samples)
  if (is.null(d)) {
    metric <- .match_metric(metric)
    d <- distance_matrix(taxa, verbose = verbose, progress = progress, metric = metric)
  }
  full <- .distinctness_matrix(d, taxa)
  method <- attr(d, "method", exact = TRUE)
  attr(full, "metric") <- if (is.null(method)) NA_character_ else method
  full
}

# MPD or MNTD of one sample (abundance vector `x` aligned with the rows of
# `full`), following picante's definitions. NA with fewer than two taxa or
# non-finite distances.
.community_statistic <- function(x, full, index, abundance_weighted) {
  present <- x > 0
  if (sum(present) < 2L) return(NA_real_)
  w <- full[present, present, drop = FALSE]
  a <- x[present]
  if (index == "mpd") {
    value <- if (abundance_weighted) {
      sum(w * outer(a, a)) / sum(outer(a, a))
    } else {
      mean(w[lower.tri(w)])
    }
  } else {
    diag(w) <- Inf
    nearest <- apply(w, 1L, min)
    value <- if (abundance_weighted) sum(nearest * a) / sum(a) else mean(nearest)
  }
  if (!is.finite(value)) NA_real_ else value
}

# Independent swap randomization (Gotelli 2000): repeatedly picks two
# samples and two taxa forming a checkerboard (each sample has one of the two
# taxa) and swaps them, keeping every sample's richness and every taxon's
# number of occurrences. Abundances move with the presences. Stops after
# `swaps` successful swaps or 100 * swaps attempts.
.independent_swap <- function(samples, swaps) {
  m <- samples
  n_row <- nrow(m)
  n_col <- ncol(m)
  done <- 0L
  attempts <- 0L
  while (done < swaps && attempts < 100L * swaps) {
    attempts <- attempts + 1L
    r <- sample.int(n_row, 2L)
    k <- sample.int(n_col, 2L)
    a <- m[r[1L], k[1L]] > 0
    b <- m[r[1L], k[2L]] > 0
    c <- m[r[2L], k[1L]] > 0
    e <- m[r[2L], k[2L]] > 0
    if (a && e && !b && !c) {
      m[r[1L], k[2L]] <- m[r[1L], k[1L]]
      m[r[1L], k[1L]] <- 0
      m[r[2L], k[1L]] <- m[r[2L], k[2L]]
      m[r[2L], k[2L]] <- 0
      done <- done + 1L
    } else if (!a && !e && b && c) {
      m[r[1L], k[1L]] <- m[r[1L], k[2L]]
      m[r[1L], k[2L]] <- 0
      m[r[2L], k[2L]] <- m[r[2L], k[1L]]
      m[r[2L], k[1L]] <- 0
      done <- done + 1L
    }
  }
  m
}
