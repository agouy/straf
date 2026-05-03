# Parity check: TS reimplementation vs original R `straf` package

Cross-validates the TypeScript stats in [`src/stats/`](../src/stats) against
the original R package output, using a small set of edge-case datasets. Not a
routine test — it requires R to be installed alongside Node — but should be
run after any change to a stats module to catch regressions against the
reference implementation.

## Layout

```
parity/
├── README.md          — this file
├── manifest.json      — list of datasets + ploidy
├── datasets/          — TSV inputs (committed; same files for both sides)
├── dump_r.R           — runs the R package on every dataset → results_r/<name>.json
├── dump_ts.ts         — runs the TS code on every dataset   → results_ts/<name>.json
├── compare.ts         — diffs the two sets of JSONs with per-section tolerances
├── results_r/         — generated, gitignored
└── results_ts/        — generated, gitignored
```

## Running

From the repo root, in this order:

```bash
# 1. R-side dump (loads the straf R package via devtools::load_all)
Rscript straf-web/parity/dump_r.R

# 2. TS-side dump
cd straf-web && npm run parity:dump

# 3. Compare
npm run parity:compare           # summary only
npm run parity:compare -- --show-expected   # full detail of expected diffs too
```

A clean run prints `All checks passed (modulo expected differences).` and
exits 0. Any unexpected mismatch exits 1 and prints the JSON path, both
values, and the absolute diff so you can drill into the source.

## What is checked

For each dataset:

| Section | What | Tolerance |
| --- | --- | --- |
| `indices` | Forensic + popgen per-locus indices: N, Nall, GD (Hexp), PIC, PM, PD, Hobs, PE, TPI | abs/rel = 1e-6 |
| `freqs` | Per-locus allele frequencies + N | abs/rel = 1e-6 |
| `perLocusFstats` | Ht, Hs, Fis, **Fst** (per locus, multi-pop diploid only) | abs/rel = 1e-6 |
| `pairwiseFst` | Weir & Cockerham θ̂ between every pair of populations | abs/rel = 1e-6 |
| `distances` | Nei, Edwards, Reynolds, Rogers, Provesti — population × population | abs/rel = 1e-6 |
| `pca` | Variance proportions + pairwise distances on first ≤5 PCs (sign-invariant) | abs/rel = 1e-4 |
| `haplotype` | Counts, frequencies, sumP², R-style diversity (haploid only) | abs/rel = 1e-9 |

## What is *not* checked

These were intentionally excluded — the user request was "everything except
the genepop ones", and a few R-only paths don't have a meaningful TS
counterpart:

- **HWE p-values** — R uses Genepop's exact MCMC; TS uses an asymptotic χ²
  approximation. Different by design.
- **LD p-values** — same: Genepop MCMC vs χ².
- **isoMDS** — R uses MASS::isoMDS; TS has its own implementation. The
  classical (cmdscale) MDS coordinates are checked indirectly via the
  population distances they're computed from.

## Datasets

10 small datasets covering common shapes and known edge cases:

| Name | Ploidy | Why |
| --- | --- | --- |
| `tiny_single_pop_diploid` | 2 | Minimal sanity (no F-stats, no MDS) |
| `two_pop_two_loci_diploid` | 2 | Smallest valid pairwise Fst |
| `multi_pop_diploid` | 2 | Common case (3 pops, 4 loci) |
| `with_missing_diploid` | 2 | `"0"` alleles scattered through the table |
| `point_alleles_diploid` | 2 | Includes 9.3 / 11.3 fractional alleles |
| `monomorphic_diploid` | 2 | One locus is fully monomorphic — division-by-zero edge |
| `tiny_haploid` | 1 | Minimal haploid |
| `multi_pop_haploid` | 1 | Multi-pop haploid |
| `haploid_missing` | 1 | Haploid with `"0"` missing markers |
| `haploid_point_alleles` | 1 | Haploid with fractional alleles |

To add a new dataset: drop a TSV in `datasets/`, add an entry in
`manifest.json`, and re-run both dumps + compare.

## Expected differences

These are documented mismatches that the comparator filters out (see the
`EXPECTED_DIFFERENCES` array in [compare.ts](compare.ts) for the live list).
Each is a place where the two implementations *legitimately* disagree —
adding a new entry should be a deliberate decision after investigating the
underlying cause.

### 1. Haploid `"0"` handling — `haploid_missing` dataset

`createGenind` in the R package converts `"0"` to `NA` *only on the
diploid branch*; the haploid branch passes `"0"` through to
`adegenet::df2genind` as a literal allele. As a result, R reports an
extra "0" allele at every locus that contains missing data, inflating N
and Nall and skewing every downstream stat.

The TS parser treats `"0"` as missing in both ploidies, which matches the
documented input format and is internally consistent. We treat the R
behavior as a pre-existing bug and accept the divergence.

### 2. Diploid partial genotypes — `with_missing_diploid` dataset

R's `createGenind` builds the diploid count matrix with
`apply(mat == allele, 1, sum)` — and `sum(c(TRUE, NA))` is `NA`. So if an
individual has *one* allele missing at a locus (e.g. `10/0`), R drops the
whole genotype. TS counts the typed allele as a single gene copy. This
shifts N, allele frequencies, and any per-locus stat that depends on
them (and therefore Fst, distances, PCA).

Both behaviours are defensible; the TS version preserves more
information. Calling this out so it can be revisited if needed — flip the
TS parser if you want hard parity here.

### 3. Haploid point-allele encoding — `haploid_point_alleles` dataset

`createGenind` for haploid data strips the dot from point alleles
(`"9.3"` → `"93"`), pads to 3 chars (`"930"`), and passes that to
`df2genind` with `ncode = 3`. The encoding is lossy: `"100"` could be
`"10"` or `"10.0"` — there's no way to recover the original dot position
from the encoded number alone. The comparator decodes what it can and
falls back to a *value-multiset* check (`freqValuesMultiset.*`) that
compares the sorted frequency distribution per locus, ignoring labels.

### 4. PCA on a fully monomorphic column — `monomorphic_diploid` dataset

`prcomp(scale = TRUE)` divides by per-column SD; a constant column has
SD = 0 and the call fails. The original R STRAF surfaces this as an
error. TS detects zero-variance columns and zeroes them, so it returns a
PCA on the remaining columns. Recorded as expected; the TS behaviour is
strictly more permissive.

## When the run *isn't* clean

If `parity:compare` reports unexpected mismatches:

1. **Read the path.** It points exactly at the field that diverged
   (e.g. `multi_pop_diploid.indices.A.L2.GD`). Open the TS source under
   `straf-web/src/stats/` and the R helper under `R/` (e.g.
   [`R/module_for_popgen.R`](../../R/module_for_popgen.R) for indices) for
   that exact field.
2. **Inspect the JSONs.** `results_r/<dataset>.json` and
   `results_ts/<dataset>.json` are pretty-printed and self-explanatory.
3. **Decide.** Either the TS implementation needs fixing, or there's a
   genuine, deliberate behavioural difference — in which case add an
   `EXPECTED_DIFFERENCES` entry in [compare.ts](compare.ts) with a clear
   reason explaining why.
