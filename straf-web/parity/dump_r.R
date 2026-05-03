#!/usr/bin/env Rscript
# Parity dump (R side).
#
# Loads the original `straf` package and runs the canonical pipeline on every
# dataset listed in parity/manifest.json, then writes one JSON file per dataset
# under parity/results_r/.  The TS-side dump produces the same JSON layout in
# parity/results_ts/ so `compare.ts` can diff them field-by-field.
#
# Skipped on purpose:
#   - HWE p-values  (Genepop MCMC; not reimplemented)
#   - LD p-values   (Genepop MCMC; not reimplemented)
#   - isoMDS        (MASS::isoMDS; we use a custom implementation)
#
# Run from the repo root:
#   Rscript straf-web/parity/dump_r.R
#
# Requires: jsonlite, plus all the packages the straf R package depends on
# (adegenet, hierfstat, pegas, etc.).
suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(straf)
  }
  library(jsonlite)
  library(adegenet)
  library(hierfstat)
})

here   <- normalizePath("straf-web/parity", winslash = "/", mustWork = TRUE)
ds_dir <- file.path(here, "datasets")
out_dir <- file.path(here, "results_r")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

manifest <- jsonlite::fromJSON(file.path(here, "manifest.json"))$datasets

# Helper: coerce a getIndicesFromGenind data frame into a tidy list of named
# numeric vectors keyed by locus name. Drops the "locus" column and renames
# "GD (Hexp)" → "GD" so it matches the TS keys.
dump_indices <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(list())
  rownames(df) <- df$locus
  df$locus <- NULL
  if ("GD (Hexp)" %in% colnames(df)) {
    colnames(df)[colnames(df) == "GD (Hexp)"] <- "GD"
  }
  loci <- rownames(df)
  out <- list()
  for (l in loci) {
    row <- df[l, , drop = FALSE]
    cols <- colnames(row)
    vals <- as.list(unlist(row))
    names(vals) <- cols
    # Cast to numeric (R sometimes returns factor / character for mixed cols).
    vals <- lapply(vals, function(v) {
      if (is.null(v) || is.na(v)) return(NA_real_)
      suppressWarnings(as.numeric(v))
    })
    out[[l]] <- vals
  }
  out
}

# Helper: getFreqFromGenind returns a character matrix where the first column
# is "Allele" and the trailing row is "N". Reshape into per-locus { allele -> freq }
# maps + per-locus N counts.
dump_freqs <- function(mat) {
  if (is.null(mat)) return(list())
  alleles <- mat[, "Allele"]
  loc_names <- setdiff(colnames(mat), "Allele")
  by_locus <- list()
  for (l in loc_names) {
    col <- mat[, l]
    n_idx <- which(alleles == "N")
    n_val <- if (length(n_idx) == 1) suppressWarnings(as.numeric(col[n_idx])) else NA_real_
    af <- list()
    for (i in seq_along(alleles)) {
      a <- alleles[i]
      if (a == "N") next
      v <- suppressWarnings(as.numeric(col[i]))
      if (!is.na(v)) af[[a]] <- v
    }
    by_locus[[l]] <- list(N = n_val, freq = af)
  }
  by_locus
}

dump_dist <- function(g, method) {
  obj <- adegenet::genind2genpop(g, quiet = TRUE)
  d <- as.matrix(adegenet::dist.genpop(obj, method = method))
  list(
    populations = rownames(d),
    matrix = unname(lapply(seq_len(nrow(d)), function(i) as.numeric(d[i, ])))
  )
}

dump_perlocus_fstats <- function(g) {
  bs <- hierfstat::basic.stats(g, diploid = TRUE, digits = 8)$perloc
  rownames(bs) <- as.character(unique(g@loc.fac))
  wc <- hierfstat::wc(g, diploid = TRUE)$per.loc$FST
  names(wc) <- as.character(unique(g@loc.fac))
  out <- list()
  for (l in rownames(bs)) {
    out[[l]] <- list(
      Ho  = unname(suppressWarnings(as.numeric(bs[l, "Ho"]))),
      Hs  = unname(suppressWarnings(as.numeric(bs[l, "Hs"]))),
      Ht  = unname(suppressWarnings(as.numeric(bs[l, "Ht"]))),
      Fis = unname(suppressWarnings(as.numeric(bs[l, "Fis"]))),
      Fst = unname(suppressWarnings(as.numeric(wc[l])))
    )
  }
  out
}

dump_pairwise_fst <- function(g) {
  m <- hierfstat::pairwise.WCfst(hierfstat::genind2hierfstat(g))
  list(
    populations = rownames(m),
    matrix = unname(lapply(seq_len(nrow(m)), function(i) as.numeric(m[i, ])))
  )
}

dump_pca <- function(g) {
  freq_tab <- adegenet::makefreq(g, missing = "mean", quiet = TRUE)
  pc <- stats::prcomp(freq_tab, scale = TRUE, center = TRUE)
  imp <- summary(pc)$importance
  list(
    nInd       = nrow(freq_tab),
    nCols      = ncol(freq_tab),
    sdev       = as.numeric(pc$sdev),
    eigenvalues= as.numeric(pc$sdev) ^ 2,
    variance   = as.numeric(imp["Proportion of Variance", ]),
    cumulative = as.numeric(imp["Cumulative Proportion", ]),
    # Pairwise distances between individuals on the first min(5, ncomp) PCs —
    # rotation/reflection invariant, lets the comparator check scores without
    # worrying about sign flips.
    pairDist   = (function() {
      k <- min(5, ncol(pc$x))
      d <- as.matrix(stats::dist(pc$x[, seq_len(k), drop = FALSE]))
      unname(lapply(seq_len(nrow(d)), function(i) as.numeric(d[i, ])))
    })(),
    individuals = rownames(freq_tab)
  )
}

dump_haplotype <- function(g) {
  # Match the original helpers_haplo_stats.R output (per-population list).
  hs <- getHaploStatsAllPop(g)
  out <- list()
  # `data.frame(h_count = a_table)` unrolls to columns `h_count.Var1` (key)
  # and `h_count.Freq` (value); pick the right column without assuming.
  pick_col <- function(hd, prefix) {
    cands <- c(paste0(prefix, ".Freq"), prefix, "Freq")
    found <- intersect(cands, colnames(hd))
    if (length(found) == 0) {
      stop(sprintf(
        "no '%s'-like column in haplotype data.frame; have: %s",
        prefix, paste(colnames(hd), collapse = ",")
      ))
    }
    found[1]
  }
  for (pop in names(hs)) {
    h <- hs[[pop]]
    if (is.null(h)) next
    hd <- h$hap_data
    hd_count <- as.numeric(hd[[pick_col(hd, "h_count")]])
    hd_freq  <- as.numeric(hd[[pick_col(hd, "h_freq")]])
    h_div_R  <- h$hap_pop_data$h_div
    sumP2    <- sum(hd_freq ^ 2)
    n        <- sum(hd_count)
    out[[pop]] <- list(
      n               = n,
      nLoci           = ncol(adegenet::genind2df(g)) - 1,
      haplotypes      = unname(lapply(seq_len(nrow(hd)), function(i) {
        list(
          id        = as.character(hd$haplotype_id[i]),
          haplotype = as.character(hd$h_value[i]),
          count     = hd_count[i],
          frequency = hd_freq[i]
        )
      })),
      sumP2           = sumP2,
      h_div_R         = unname(as.numeric(h_div_R)),
      diversity_TSlike = if (n > 1) (n / (n - 1)) * (1 - sumP2) else NA_real_,
      d_mat = unname(lapply(seq_len(nrow(h$d_mat)),
                            function(i) as.numeric(h$d_mat[i, ])))
    )
  }
  out
}

results <- list()
for (i in seq_len(nrow(manifest))) {
  ds <- manifest[i, ]
  name   <- ds$name
  ploidy <- ds$ploidy
  fpath  <- file.path(ds_dir, paste0(name, ".txt"))
  cat(sprintf("[R] %s (ploidy=%d) ...\n", name, ploidy))

  g <- createGenind(list(datapath = fpath), ploidy = ploidy)

  out <- list(
    dataset = name,
    ploidy  = ploidy,
    nInd    = nrow(g@tab),
    pops    = sort(unique(as.character(pop(g))))
  )

  # Per-population indices + freqs.
  ind_all <- getIndicesAllPop(g, ploidy = ploidy)
  freq_all <- getFreqAllPop(g)
  out$indices <- list()
  out$freqs   <- list()
  for (p in names(ind_all)) {
    out$indices[[p]] <- dump_indices(ind_all[[p]])
  }
  for (p in names(freq_all)) {
    out$freqs[[p]] <- dump_freqs(freq_all[[p]])
  }

  # Multi-pop diploid extras: F-stats, pairwise Fst, distances.
  npop <- length(unique(g@pop))
  if (ploidy == 2 && npop > 1 && length(locNames(g)) > 1) {
    out$perLocusFstats <- tryCatch(dump_perlocus_fstats(g), error = function(e) NULL)
    out$pairwiseFst    <- tryCatch(dump_pairwise_fst(g),    error = function(e) NULL)
    out$distances <- list()
    methods <- list(nei = 1, edwards = 2, reynolds = 3, rogers = 4, provesti = 5)
    for (mname in names(methods)) {
      out$distances[[mname]] <- tryCatch(dump_dist(g, methods[[mname]]),
                                         error = function(e) NULL)
    }
  }

  # PCA: needs ≥ 2 individuals.
  if (nrow(g@tab) >= 2) {
    out$pca <- tryCatch(dump_pca(g), error = function(e) NULL)
  }

  # Haploid: haplotype stats.
  if (ploidy == 1) {
    out$haplotype <- tryCatch(dump_haplotype(g), error = function(e) NULL)
  }

  # File conversions: capture the raw text of each output. Some conversions
  # don't apply to all ploidies (R's straf2arlequin pairs columns assuming
  # diploid; straf2familias has the same assumption).
  out$conversions <- list()

  out$conversions$genepop <- tryCatch({
    tmp <- tempfile(); on.exit(unlink(tmp), add = TRUE)
    straf2genepop(fpath, tmp, ploidy = ploidy)
    paste(readLines(tmp, warn = FALSE), collapse = "\n")
  }, error = function(e) NULL)

  if (ploidy == 2) {
    out$conversions$arlequin <- tryCatch({
      tmp <- tempfile(); on.exit(unlink(tmp), add = TRUE)
      straf2arlequin(fpath, tmp)
      paste(readLines(tmp, warn = FALSE), collapse = "\n")
    }, error = function(e) NULL)

    out$conversions$familias <- list()
    pops_for_freq <- c("all", out$pops)
    for (p in pops_for_freq) {
      out$conversions$familias[[p]] <- tryCatch({
        tmp <- tempfile(); on.exit(unlink(tmp), add = TRUE)
        straf2familias(fpath, tmp, pop = p)
        paste(readLines(tmp, warn = FALSE), collapse = "\n")
      }, error = function(e) NULL)
    }
  }

  # Euroformix / STRmix / LRmix: these pull from getFreqAllPop, which works
  # for both ploidies. We mirror the UI's write.table calls exactly.
  freqs_for_csv <- tryCatch(getFreqAllPop(g), error = function(e) NULL)
  if (!is.null(freqs_for_csv)) {
    write_csv_variant <- function(matr, drop_n_row, na_repr) {
      tmp <- tempfile(); on.exit(unlink(tmp), add = TRUE)
      m <- if (drop_n_row) matr[-nrow(matr), , drop = FALSE] else matr
      utils::write.table(
        m, tmp, sep = ",",
        na = na_repr, row.names = FALSE, quote = FALSE
      )
      paste(readLines(tmp, warn = FALSE), collapse = "\n")
    }
    out$conversions$euroformix <- list()
    out$conversions$strmix <- list()
    out$conversions$lrmix <- list()
    for (p in names(freqs_for_csv)) {
      matr <- freqs_for_csv[[p]]
      out$conversions$euroformix[[p]] <- tryCatch(
        write_csv_variant(matr, drop_n_row = TRUE, na_repr = ""),
        error = function(e) NULL
      )
      out$conversions$strmix[[p]] <- tryCatch(
        write_csv_variant(matr, drop_n_row = FALSE, na_repr = "0"),
        error = function(e) NULL
      )
      out$conversions$lrmix[[p]] <- tryCatch(
        write_csv_variant(matr, drop_n_row = TRUE, na_repr = ""),
        error = function(e) NULL
      )
    }
  }

  jsonlite::write_json(
    out,
    file.path(out_dir, paste0(name, ".json")),
    auto_unbox = TRUE,
    pretty     = TRUE,
    digits     = NA,
    null       = "null",
    na         = "null"
  )
  results[[name]] <- "ok"
}

cat(sprintf("[R] wrote %d JSON files to %s\n", length(results), out_dir))
