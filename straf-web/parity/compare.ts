/**
 * Parity comparator.
 *
 * Walks the JSON files written by `dump_r.R` and `dump_ts.ts`, diffing
 * numeric fields with section-specific tolerances. Reports each mismatch with
 * a path like `multi_pop_diploid.indices.A.L2.GD` so it's obvious where to
 * look. Exits with a non-zero status if any check fails.
 *
 * Run from straf-web/:
 *   node --import tsx parity/compare.ts
 */

import { readFileSync, readdirSync } from "node:fs";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const rDir = join(here, "results_r");
const tsDir = join(here, "results_ts");

interface Tolerance {
  /** Absolute tolerance: fail if |a − b| > abs. */
  abs: number;
  /** Relative tolerance: fail if |a − b| > rel · max(|a|, |b|). */
  rel: number;
}

const DEFAULT_TOL: Tolerance = { abs: 1e-6, rel: 1e-6 };
const PCA_TOL: Tolerance = { abs: 1e-4, rel: 1e-4 };
const HAPLO_DIV_TOL: Tolerance = { abs: 1e-9, rel: 1e-9 };

interface Mismatch {
  dataset: string;
  path: string;
  r: unknown;
  ts: unknown;
  diff?: number;
  reason: string;
  expected?: string; // populated by classifyMismatches if it matches a known difference
}

const mismatches: Mismatch[] = [];

function recordMismatch(m: Mismatch): void {
  mismatches.push(m);
}

/**
 * Known, intentional or pre-existing differences between the R package output
 * and the TS reimplementation. Mismatches that match any entry below are
 * filtered out of the failure count and reported in a separate "expected"
 * section. Add a new entry only after explicitly deciding the difference is
 * acceptable — every entry here is a place where the two implementations
 * legitimately disagree.
 */
interface ExpectedDifference {
  /** Dataset name pattern (substring match). */
  dataset?: RegExp;
  /** JSON-path pattern (regex match). */
  path: RegExp;
  /** Why this difference is expected. */
  reason: string;
}

const EXPECTED_DIFFERENCES: ExpectedDifference[] = [
  // ── Haploid: R's createGenind for ploidy=1 doesn't convert "0" to NA
  // (only the diploid branch does). Result: "0" gets treated as a real
  // allele in R's haploid path. The TS parser treats it as missing in both
  // ploidies (consistent + matches the documented input format).
  {
    dataset: /^haploid_missing$/,
    path: /^(indices|freqs)\./,
    reason: "R haploid path treats \"0\" as a real allele (not converted to NA). TS treats it as missing.",
  },
  {
    dataset: /^haploid_missing$/,
    path: /^haplotype\./,
    reason: "Same root cause: R haploid \"0\" handling.",
  },
  {
    dataset: /^haploid_missing$/,
    path: /^pca\./,
    reason: "Same root cause: R haploid \"0\" handling shifts the dosage matrix.",
  },
  {
    dataset: /^haploid_missing$/,
    path: /^conversions\.(genepop|euroformix|strmix|lrmix)\./,
    reason: "Same root cause: R haploid \"0\" handling propagates into all freq-based conversion outputs.",
  },

  // ── Diploid with missing alleles: R drops the entire genotype at a locus
  // when *either* allele is missing (apply+sum on NA → NA). TS counts the
  // typed allele (so a 10/0 contributes 1 to N at that locus). This affects
  // N, allele frequencies and any per-locus stat that depends on them.
  {
    dataset: /^with_missing_diploid$/,
    path: /^(indices|freqs)\./,
    reason: "R drops partial diploid genotypes (one allele missing → both treated as missing). TS counts the typed allele.",
  },
  {
    dataset: /^with_missing_diploid$/,
    path: /^perLocusFstats\.[^.]+\.Fst$/,
    reason: "Knock-on effect of R's partial-genotype dropping.",
  },
  {
    dataset: /^with_missing_diploid$/,
    path: /^pairwiseFst\./,
    reason: "Knock-on effect of R's partial-genotype dropping.",
  },
  {
    dataset: /^with_missing_diploid$/,
    path: /^distances\./,
    reason: "Knock-on effect of R's partial-genotype dropping.",
  },
  {
    dataset: /^with_missing_diploid$/,
    path: /^pca\./,
    reason: "Knock-on effect of R's partial-genotype dropping.",
  },
  {
    dataset: /^with_missing_diploid$/,
    path: /^conversions\.familias\./,
    reason: "R's straf2familias reads the file directly and counts \"0\" as a real allele (raw table() on strings, no genind round-trip). TS doesn't.",
  },
  {
    dataset: /^with_missing_diploid$/,
    path: /^conversions\.(euroformix|strmix|lrmix)\./,
    reason: "Knock-on of R dropping partial diploid genotypes when building the genind that feeds the freq table.",
  },

  // ── PCA on a dataset with a monomorphic locus: prcomp(scale=TRUE) fails
  // on a constant column (sd = 0 → division by zero). Our TS PCA handles
  // it by zeroing the column. So R produces no PCA output for this dataset.
  {
    dataset: /^monomorphic_diploid$/,
    path: /^pca/,
    reason: "R's prcomp(scale=TRUE) fails on a constant (monomorphic) column. TS handles it.",
  },

  // ── R encodes haploid alleles as 3-char strings ("9.3" → "930"); leading
  // zeros and the dot position are lost when as.numeric() coerces them
  // back to numbers for the freq table key. The encoding is lossy ("100"
  // could be "10" or "10.0"), so labels can't be matched up reliably. We
  // fall back to a value-multiset check (see freqValuesMatch) which
  // verifies the *distribution* of frequencies per locus matches.
  {
    dataset: /^haploid_point_alleles$/,
    path: /^freqs\./,
    reason: "R haploid encodes alleles as 3-char strings; dot position is lost. Value-multiset check still runs.",
  },
  {
    dataset: /^haploid_point_alleles$/,
    path: /^haplotype\.[^.]+\.haplotypes\[\d+\]\.haplotype$/,
    reason: "Haplotype string concatenates encoded allele codes — same encoding-lossiness as freqs.",
  },
  // Same root cause: per-allele label mismatches in haploid CSV exports.
  // The compareCsvOutputs function then runs a sorted-multiset value check
  // (`values_sorted`) which is the source of truth for these.
  {
    dataset: /^haploid_point_alleles$/,
    path: /^conversions\.(euroformix|strmix|lrmix)\..*\.[^.]*$/,
    reason: "Same root cause as above; per-CSV multiset value check confirms the underlying distribution.",
  },

  // ── Arlequin output: R's straf2arlequin has a known encoding bug. When no
  // allele in a locus contains a dot, alleles are emitted as raw strings
  // ("10"); when any allele has a dot, the encoding pads inconsistently to
  // 2-3 chars. TS produces the spec-compliant 3-digit form everywhere.
  // Recorded as expected — TS is correct here.
  {
    path: /^conversions\.arlequin\./,
    reason: "R's straf2arlequin produces variable-length allele codes (encoding bug). TS uses spec-compliant 3-digit form.",
  },
];

function classifyMismatches(): void {
  for (const m of mismatches) {
    for (const e of EXPECTED_DIFFERENCES) {
      if (e.dataset && !e.dataset.test(m.dataset)) continue;
      if (!e.path.test(m.path)) continue;
      m.expected = e.reason;
      break;
    }
  }
}

/**
 * Inverse of R's createGenind haploid encoding: a 3-char numeric code where
 * the leading char is a "0"-pad, the middle/last is the integer, and a
 * fractional digit is the trailing char (with "0" meaning "no fraction").
 *
 *   "090" → "9"     (no fraction, padded)
 *   "100" → "10"    (no fraction)
 *   "930" → "9.3"   (fraction)
 *
 * Used for the haploid_point_alleles dataset where R-side allele keys won't
 * match TS keys. The encoding is described inline in createGenind:
 *   gsub("[.]", "", x);
 *   nchar(x)==1 → paste0("0", x, "0");
 *   nchar(x)==2 → paste0(x, "0");
 *   else require nchar(x)==3.
 */
function decodeHaploidAllele(code: string): string {
  if (!/^[0-9]{3}$/.test(code)) return code;
  const intPart = String(parseInt(code.slice(0, 2), 10));
  const frac = code.charAt(2);
  return frac === "0" ? intPart : `${intPart}.${frac}`;
}

function isHaploidPointDataset(dataset: string): boolean {
  return /haploid.*point|point.*haploid/i.test(dataset);
}

/**
 * Per-locus, per-population: compare the *sorted multiset* of allele
 * frequencies (and N) between R and TS, ignoring the allele *labels*. Used
 * for haploid_point_alleles where R's encoding loses dot positions but the
 * underlying frequency distribution should be identical.
 */
function compareFreqValuesMultiset(
  dataset: string,
  rF: Record<string, Record<string, { N: number; freq: Record<string, number> }>>,
  tsF: Record<string, Record<string, { N: number; freq: Record<string, number> }>>,
): void {
  for (const pop of new Set([...Object.keys(rF), ...Object.keys(tsF)])) {
    const rPop = rF[pop];
    const tsPop = tsF[pop];
    if (!rPop || !tsPop) continue;
    for (const loc of new Set([...Object.keys(rPop), ...Object.keys(tsPop)])) {
      const rLoc = rPop[loc];
      const tsLoc = tsPop[loc];
      if (!rLoc || !tsLoc) continue;
      compareNumbers(dataset, `freqValuesMultiset.${pop}.${loc}.N`, rLoc.N, tsLoc.N, DEFAULT_TOL);
      const rVals = Object.values(rLoc.freq).sort((a, b) => a - b);
      const tsVals = Object.values(tsLoc.freq).sort((a, b) => a - b);
      compareNumberArray(dataset, `freqValuesMultiset.${pop}.${loc}`, rVals, tsVals, DEFAULT_TOL);
    }
  }
}

function remapAlleles(
  src: Record<string, Record<string, { N: number; freq: Record<string, number> }>>,
  decode: (s: string) => string,
): Record<string, Record<string, { N: number; freq: Record<string, number> }>> {
  const out: typeof src = {};
  for (const pop of Object.keys(src)) {
    out[pop] = {};
    for (const loc of Object.keys(src[pop]!)) {
      const inFreq = src[pop]![loc]!.freq;
      const newFreq: Record<string, number> = {};
      for (const a of Object.keys(inFreq)) newFreq[decode(a)] = inFreq[a]!;
      out[pop]![loc] = { N: src[pop]![loc]!.N, freq: newFreq };
    }
  }
  return out;
}

function nearlyEqual(a: number, b: number, tol: Tolerance): boolean {
  if (Number.isNaN(a) && Number.isNaN(b)) return true;
  if (!Number.isFinite(a) || !Number.isFinite(b)) return a === b;
  const diff = Math.abs(a - b);
  if (diff <= tol.abs) return true;
  return diff <= tol.rel * Math.max(Math.abs(a), Math.abs(b));
}

function compareNumbers(
  dataset: string,
  path: string,
  r: unknown,
  ts: unknown,
  tol: Tolerance,
): void {
  // Treat null/undefined as NaN — both should agree on missingness.
  const rNull = r === null || r === undefined;
  const tsNull = ts === null || ts === undefined;
  if (rNull && tsNull) return;
  if (rNull || tsNull) {
    recordMismatch({ dataset, path, r, ts, reason: "one side is null/undefined" });
    return;
  }
  const rn = typeof r === "number" ? r : Number(r);
  const tn = typeof ts === "number" ? ts : Number(ts);
  if (Number.isNaN(rn) && Number.isNaN(tn)) return;
  if (!Number.isFinite(rn) && !Number.isFinite(tn)) {
    if (rn === tn || (Number.isNaN(rn) && Number.isNaN(tn))) return;
    recordMismatch({ dataset, path, r, ts, reason: "non-finite mismatch" });
    return;
  }
  if (!nearlyEqual(rn, tn, tol)) {
    recordMismatch({
      dataset,
      path,
      r: rn,
      ts: tn,
      diff: Math.abs(rn - tn),
      reason: `out of tolerance (abs=${tol.abs}, rel=${tol.rel})`,
    });
  }
}

function compareNumberArray(
  dataset: string,
  path: string,
  r: unknown,
  ts: unknown,
  tol: Tolerance,
): void {
  if (!Array.isArray(r) || !Array.isArray(ts)) {
    recordMismatch({ dataset, path, r, ts, reason: "expected arrays" });
    return;
  }
  if (r.length !== ts.length) {
    recordMismatch({ dataset, path, r: r.length, ts: ts.length, reason: "array length differs" });
    return;
  }
  for (let i = 0; i < r.length; i++) {
    compareNumbers(dataset, `${path}[${i}]`, r[i], ts[i], tol);
  }
}

function compareNumberMatrix(
  dataset: string,
  path: string,
  r: unknown,
  ts: unknown,
  tol: Tolerance,
): void {
  if (!Array.isArray(r) || !Array.isArray(ts)) {
    recordMismatch({ dataset, path, r, ts, reason: "expected matrices" });
    return;
  }
  if (r.length !== ts.length) {
    recordMismatch({ dataset, path, r: r.length, ts: ts.length, reason: "row count differs" });
    return;
  }
  for (let i = 0; i < r.length; i++) {
    compareNumberArray(dataset, `${path}[${i}]`, r[i], ts[i], tol);
  }
}

/**
 * Reorder the rows/columns of `tsMatrix` (with row labels `tsLabels`) so that
 * they match `rLabels`. Returns null if the label sets differ.
 */
function reorderSquareMatrix(
  tsLabels: string[],
  tsMatrix: number[][],
  rLabels: string[],
): number[][] | null {
  if (tsLabels.length !== rLabels.length) return null;
  const idx = new Map(tsLabels.map((l, i) => [l, i] as const));
  for (const l of rLabels) if (!idx.has(l)) return null;
  return rLabels.map((rl) => {
    const ri = idx.get(rl)!;
    return rLabels.map((cl) => tsMatrix[ri]![idx.get(cl)!]!);
  });
}

function comparePerPopIndices(
  dataset: string,
  rIdx: Record<string, Record<string, Record<string, unknown>>>,
  tsIdx: Record<string, Record<string, Record<string, unknown>>>,
): void {
  const pops = new Set<string>([...Object.keys(rIdx), ...Object.keys(tsIdx)]);
  for (const pop of pops) {
    const rPop = rIdx[pop];
    const tsPop = tsIdx[pop];
    if (!rPop || !tsPop) {
      recordMismatch({
        dataset,
        path: `indices.${pop}`,
        r: !!rPop,
        ts: !!tsPop,
        reason: "population missing on one side",
      });
      continue;
    }
    const loci = new Set<string>([...Object.keys(rPop), ...Object.keys(tsPop)]);
    for (const loc of loci) {
      const rRow = rPop[loc];
      const tsRow = tsPop[loc];
      if (!rRow || !tsRow) {
        recordMismatch({
          dataset,
          path: `indices.${pop}.${loc}`,
          r: !!rRow,
          ts: !!tsRow,
          reason: "locus missing on one side",
        });
        continue;
      }
      // Only compare keys present on both sides — TS exposes Hobs/PE/TPI for
      // diploid, R adds Ht/Fis/Fst when ≥ 2 pops; we compare those separately.
      const sharedKeys = ["N", "Nall", "GD", "PIC", "PM", "PD", "Hobs", "PE", "TPI"];
      for (const k of sharedKeys) {
        if (!(k in rRow) || !(k in tsRow)) continue;
        compareNumbers(dataset, `indices.${pop}.${loc}.${k}`, rRow[k], tsRow[k], DEFAULT_TOL);
      }
    }
  }
}

function comparePerPopFreqs(
  dataset: string,
  rF: Record<string, Record<string, { N: number; freq: Record<string, number> }>>,
  tsF: Record<string, Record<string, { N: number; freq: Record<string, number> }>>,
): void {
  // R's haploid path encodes alleles as 3-char strings ("9.3" → "930"); decode
  // so allele keys can be matched up against the TS side. The encoding is
  // lossy (`100` could be `10` or `10.0`), so for that dataset we *also* run a
  // value-multiset check below as the source of truth.
  if (isHaploidPointDataset(dataset)) {
    rF = remapAlleles(rF, decodeHaploidAllele);
    compareFreqValuesMultiset(dataset, rF, tsF);
  }
  const pops = new Set<string>([...Object.keys(rF), ...Object.keys(tsF)]);
  for (const pop of pops) {
    const rPop = rF[pop];
    const tsPop = tsF[pop];
    if (!rPop || !tsPop) {
      recordMismatch({
        dataset,
        path: `freqs.${pop}`,
        r: !!rPop,
        ts: !!tsPop,
        reason: "population missing on one side",
      });
      continue;
    }
    const loci = new Set<string>([...Object.keys(rPop), ...Object.keys(tsPop)]);
    for (const loc of loci) {
      const rLoc = rPop[loc];
      const tsLoc = tsPop[loc];
      if (!rLoc || !tsLoc) {
        recordMismatch({
          dataset,
          path: `freqs.${pop}.${loc}`,
          r: !!rLoc,
          ts: !!tsLoc,
          reason: "locus missing on one side",
        });
        continue;
      }
      compareNumbers(dataset, `freqs.${pop}.${loc}.N`, rLoc.N, tsLoc.N, DEFAULT_TOL);
      const alleles = new Set<string>([...Object.keys(rLoc.freq), ...Object.keys(tsLoc.freq)]);
      for (const a of alleles) {
        const rv = rLoc.freq[a];
        const tv = tsLoc.freq[a];
        // Allow allele to be absent on one side iff the frequency is 0 on
        // both sides (R sometimes omits zero entries; we sometimes do too).
        if (rv === undefined && (tv === undefined || tv === 0)) continue;
        if (tv === undefined && (rv === undefined || rv === 0)) continue;
        compareNumbers(dataset, `freqs.${pop}.${loc}.freq.${a}`, rv, tv, DEFAULT_TOL);
      }
    }
  }
}

function comparePerLocusFstats(
  dataset: string,
  r: Record<string, Record<string, number>>,
  ts: Record<string, Record<string, number>>,
): void {
  const loci = new Set<string>([...Object.keys(r), ...Object.keys(ts)]);
  for (const loc of loci) {
    const rL = r[loc];
    const tsL = ts[loc];
    if (!rL || !tsL) {
      recordMismatch({
        dataset,
        path: `perLocusFstats.${loc}`,
        r: !!rL,
        ts: !!tsL,
        reason: "locus missing on one side",
      });
      continue;
    }
    for (const k of ["Ho", "Ht", "Hs", "Fis", "Fst"]) {
      if (!(k in rL) || !(k in tsL)) continue;
      compareNumbers(dataset, `perLocusFstats.${loc}.${k}`, rL[k], tsL[k], DEFAULT_TOL);
    }
  }
}

function comparePairwiseFst(
  dataset: string,
  r: { populations: string[]; matrix: (number | null)[][] },
  ts: { populations: string[]; matrix: (number | null)[][] },
): void {
  if (!r || !ts) {
    if (r || ts) {
      recordMismatch({ dataset, path: "pairwiseFst", r: !!r, ts: !!ts, reason: "one side missing" });
    }
    return;
  }
  // hierfstat::pairwise.WCfst returns NA on the diagonal; we put 0. Replace
  // NaNs/nulls with 0 to keep the matrix comparison clean.
  const norm = (m: (number | null)[][]) =>
    m.map((row) => row.map((v) => (v === null || (typeof v === "number" && !Number.isFinite(v)) ? 0 : v)));
  const rN = norm(r.matrix);
  const tsReordered = reorderSquareMatrix(ts.populations, ts.matrix as number[][], r.populations);
  if (!tsReordered) {
    recordMismatch({
      dataset,
      path: "pairwiseFst.populations",
      r: r.populations,
      ts: ts.populations,
      reason: "population labels differ",
    });
    return;
  }
  compareNumberMatrix(dataset, "pairwiseFst.matrix", rN, norm(tsReordered as (number | null)[][]), DEFAULT_TOL);
}

function compareDistances(
  dataset: string,
  r: Record<string, { populations: string[]; matrix: number[][] }>,
  ts: Record<string, { populations: string[]; matrix: number[][] }>,
): void {
  for (const m of ["nei", "edwards", "reynolds", "rogers", "provesti"] as const) {
    const rD = r[m];
    const tD = ts[m];
    if (!rD || !tD) {
      if (rD || tD) {
        recordMismatch({
          dataset,
          path: `distances.${m}`,
          r: !!rD,
          ts: !!tD,
          reason: "one side missing",
        });
      }
      continue;
    }
    const tsReordered = reorderSquareMatrix(tD.populations, tD.matrix, rD.populations);
    if (!tsReordered) {
      recordMismatch({
        dataset,
        path: `distances.${m}.populations`,
        r: rD.populations,
        ts: tD.populations,
        reason: "population labels differ",
      });
      continue;
    }
    compareNumberMatrix(dataset, `distances.${m}.matrix`, rD.matrix, tsReordered, DEFAULT_TOL);
  }
}

function comparePca(
  dataset: string,
  r: { variance: number[]; eigenvalues: number[]; pairDist: number[][]; individuals: string[] },
  ts: { variance: number[]; eigenvalues: number[]; pairDist: number[][]; individuals: string[] },
): void {
  if (!r || !ts) {
    if (r || ts) {
      recordMismatch({ dataset, path: "pca", r: !!r, ts: !!ts, reason: "one side missing" });
    }
    return;
  }
  // Variance proportions: should match exactly (rotation/sign invariant).
  // Drop trailing zero components (numerical noise) before comparing.
  const trimZero = (a: number[]) => {
    let last = a.length - 1;
    while (last >= 0 && Math.abs(a[last]!) < 1e-10) last--;
    return a.slice(0, last + 1);
  };
  const rv = trimZero(r.variance);
  const tv = trimZero(ts.variance);
  const k = Math.min(rv.length, tv.length);
  for (let i = 0; i < k; i++) {
    compareNumbers(dataset, `pca.variance[${i}]`, rv[i], tv[i], PCA_TOL);
  }
  // Reorder TS pairwise-distance matrix to match R's individual order.
  const tsReordered = reorderSquareMatrix(ts.individuals, ts.pairDist, r.individuals);
  if (!tsReordered) {
    recordMismatch({
      dataset,
      path: "pca.individuals",
      r: r.individuals.length,
      ts: ts.individuals.length,
      reason: "individual order/labels differ",
    });
    return;
  }
  compareNumberMatrix(dataset, "pca.pairDist", r.pairDist, tsReordered, PCA_TOL);
}

function compareHaplotype(
  dataset: string,
  r: Record<string, Record<string, unknown>>,
  ts: Record<string, Record<string, unknown>>,
): void {
  const pops = new Set<string>([...Object.keys(r), ...Object.keys(ts)]);
  for (const pop of pops) {
    const rPop = r[pop] as undefined | {
      n: number;
      sumP2: number;
      h_div_R: number;
      haplotypes: { count: number; frequency: number; haplotype: string }[];
    };
    const tsPop = ts[pop] as typeof rPop;
    if (!rPop || !tsPop) {
      recordMismatch({
        dataset,
        path: `haplotype.${pop}`,
        r: !!rPop,
        ts: !!tsPop,
        reason: "population missing on one side",
      });
      continue;
    }
    compareNumbers(dataset, `haplotype.${pop}.n`, rPop.n, tsPop.n, DEFAULT_TOL);
    compareNumbers(dataset, `haplotype.${pop}.sumP2`, rPop.sumP2, tsPop.sumP2, HAPLO_DIV_TOL);
    compareNumbers(dataset, `haplotype.${pop}.h_div_R`, rPop.h_div_R, tsPop.h_div_R, HAPLO_DIV_TOL);
    // Compare per-haplotype counts/frequencies after sorting by haplotype value
    // (the two sides may emit them in slightly different orders).
    const sortKey = (h: { haplotype: string }) => h.haplotype;
    const rHaps = [...rPop.haplotypes].sort((a, b) => sortKey(a).localeCompare(sortKey(b)));
    const tsHaps = [...tsPop.haplotypes].sort((a, b) => sortKey(a).localeCompare(sortKey(b)));
    if (rHaps.length !== tsHaps.length) {
      recordMismatch({
        dataset,
        path: `haplotype.${pop}.haplotypes.length`,
        r: rHaps.length,
        ts: tsHaps.length,
        reason: "distinct haplotype count differs",
      });
      continue;
    }
    for (let i = 0; i < rHaps.length; i++) {
      const rh = rHaps[i]!;
      const th = tsHaps[i]!;
      if (rh.haplotype !== th.haplotype) {
        recordMismatch({
          dataset,
          path: `haplotype.${pop}.haplotypes[${i}].haplotype`,
          r: rh.haplotype,
          ts: th.haplotype,
          reason: "haplotype value differs",
        });
        continue;
      }
      compareNumbers(dataset, `haplotype.${pop}.haplotypes[${i}].count`, rh.count, th.count, DEFAULT_TOL);
      compareNumbers(dataset, `haplotype.${pop}.haplotypes[${i}].frequency`, rh.frequency, th.frequency, DEFAULT_TOL);
    }
  }
}

// --------- Conversion (file format) parsers ------------------------------
//
// Both implementations produce text output. We can't byte-compare because:
//   - R appends " \n" between lines (paste(x, "\n", collapse=""))
//   - R Familias uses table()-based lex order; TS sorts numerically
//   - R Arlequin has a known encoding bug (variable-length allele codes when
//     point alleles are present, raw strings when not)
// Instead, we parse each output back to a canonical structure and compare
// structurally.

interface CanonGenepop {
  loci: string[];
  pops: { individuals: { id: string; encoded: string[] }[] }[];
}

function parseGenepop(text: string): CanonGenepop {
  const lines = text
    .split(/\r?\n/)
    .map((l) => l.replace(/[\t ]+$/, ""))
    .filter((l) => l.length > 0);
  const loci: string[] = [];
  const pops: CanonGenepop["pops"] = [];
  let i = 1; // skip header
  while (i < lines.length && lines[i] !== "Pop") loci.push(lines[i++]!);
  let cur: CanonGenepop["pops"][number] | null = null;
  for (; i < lines.length; i++) {
    const l = lines[i]!;
    if (l === "Pop") {
      cur = { individuals: [] };
      pops.push(cur);
      continue;
    }
    if (!cur) continue;
    // "id\t,\t<g1>\t<g2>..."
    const parts = l.split(/\t/);
    const id = parts[0]!;
    let started = false;
    const encoded: string[] = [];
    for (let k = 1; k < parts.length; k++) {
      if (!started) {
        if (parts[k] === ",") started = true;
        continue;
      }
      if (parts[k]!.length > 0) encoded.push(parts[k]!);
    }
    cur.individuals.push({ id, encoded });
  }
  return { loci, pops };
}

function compareGenepopOutputs(dataset: string, rText: string, tsText: string): void {
  const r = parseGenepop(rText);
  const ts = parseGenepop(tsText);
  if (r.loci.length !== ts.loci.length) {
    recordMismatch({
      dataset,
      path: "conversions.genepop.loci.length",
      r: r.loci.length,
      ts: ts.loci.length,
      reason: "differs",
    });
    return;
  }
  for (let l = 0; l < r.loci.length; l++) {
    if (r.loci[l] !== ts.loci[l]) {
      recordMismatch({
        dataset,
        path: `conversions.genepop.loci[${l}]`,
        r: r.loci[l],
        ts: ts.loci[l],
        reason: "locus name differs",
      });
    }
  }
  if (r.pops.length !== ts.pops.length) {
    recordMismatch({
      dataset,
      path: "conversions.genepop.pops.length",
      r: r.pops.length,
      ts: ts.pops.length,
      reason: "differs",
    });
    return;
  }
  for (let p = 0; p < r.pops.length; p++) {
    const rPop = r.pops[p]!;
    const tsPop = ts.pops[p]!;
    const tsById = new Map(tsPop.individuals.map((ind) => [ind.id, ind] as const));
    for (const rInd of rPop.individuals) {
      const tsInd = tsById.get(rInd.id);
      if (!tsInd) {
        recordMismatch({
          dataset,
          path: `conversions.genepop.pops[${p}].${rInd.id}`,
          r: "present",
          ts: "missing",
          reason: "individual missing in TS output",
        });
        continue;
      }
      if (rInd.encoded.length !== tsInd.encoded.length) {
        recordMismatch({
          dataset,
          path: `conversions.genepop.pops[${p}].${rInd.id}.encoded.length`,
          r: rInd.encoded.length,
          ts: tsInd.encoded.length,
          reason: "encoded genotype count differs",
        });
        continue;
      }
      for (let l = 0; l < rInd.encoded.length; l++) {
        if (rInd.encoded[l] !== tsInd.encoded[l]) {
          recordMismatch({
            dataset,
            path: `conversions.genepop.pops[${p}].${rInd.id}.${r.loci[l]}`,
            r: rInd.encoded[l],
            ts: tsInd.encoded[l],
            reason: "encoded genotype differs",
          });
        }
      }
    }
  }
}

/** Familias: blocks separated by blank line. Each block: locus name + lines of "<allele>\t<freq>". */
function parseFamilias(text: string): Record<string, Record<string, number>> {
  const blocks = text
    .split(/\n[ \t]*\n/)
    .map((b) => b.trim())
    .filter((b) => b.length > 0);
  const out: Record<string, Record<string, number>> = {};
  for (const block of blocks) {
    const lines = block.split(/\r?\n/).map((l) => l.replace(/[\t ]+$/, ""));
    if (lines.length === 0) continue;
    const locus = lines[0]!;
    const freq: Record<string, number> = {};
    for (let i = 1; i < lines.length; i++) {
      const parts = lines[i]!.split("\t");
      if (parts.length < 2) continue;
      const a = parts[0]!;
      const f = Number(parts[1]!);
      if (Number.isFinite(f)) freq[a] = f;
    }
    out[locus] = freq;
  }
  return out;
}

function compareFamiliasOutputs(
  dataset: string,
  pathPrefix: string,
  rText: string,
  tsText: string,
): void {
  const r = parseFamilias(rText);
  const ts = parseFamilias(tsText);
  const loci = new Set([...Object.keys(r), ...Object.keys(ts)]);
  for (const loc of loci) {
    const rL = r[loc];
    const tsL = ts[loc];
    if (!rL || !tsL) {
      recordMismatch({
        dataset,
        path: `${pathPrefix}.${loc}`,
        r: !!rL,
        ts: !!tsL,
        reason: "locus block missing on one side",
      });
      continue;
    }
    const alleles = new Set([...Object.keys(rL), ...Object.keys(tsL)]);
    for (const a of alleles) {
      compareNumbers(dataset, `${pathPrefix}.${loc}.${a}`, rL[a], tsL[a], DEFAULT_TOL);
    }
  }
}

/**
 * Allele-frequency CSV (Euroformix / STRmix / LRmix): header row then
 * `<allele>,<freq_L1>,<freq_L2>,...` rows, optionally with a trailing N row.
 * Returns canonical structure: { locus → { allele → freq }, N: number[] }.
 *
 * Variants encode missing differently ("" for Euroformix/LRmix, "0" for
 * STRmix data rows). We treat both "" and "0" as null, since alleles with
 * literal frequency 0 are never emitted by either side.
 */
interface CanonCsv {
  loci: string[];
  freqs: Record<string, Record<string, number>>; // locus → allele → freq
  N: Record<string, number | null>; // locus → N (or null if N row absent)
}

function parseCsv(text: string): CanonCsv {
  const lines = text
    .split(/\r?\n/)
    .map((l) => l.replace(/[\t ]+$/, ""))
    .filter((l) => l.length > 0);
  if (lines.length === 0) return { loci: [], freqs: {}, N: {} };
  const headers = lines[0]!.split(",");
  const loci = headers.slice(1);
  const freqs: Record<string, Record<string, number>> = {};
  const N: Record<string, number | null> = {};
  for (const l of loci) {
    freqs[l] = {};
    N[l] = null;
  }
  for (let i = 1; i < lines.length; i++) {
    const cells = lines[i]!.split(",");
    const allele = cells[0]!;
    if (allele === "N") {
      for (let l = 0; l < loci.length; l++) {
        const v = cells[l + 1] ?? "";
        const n = Number(v);
        N[loci[l]!] = Number.isFinite(n) ? n : null;
      }
      continue;
    }
    for (let l = 0; l < loci.length; l++) {
      const v = cells[l + 1] ?? "";
      if (v === "" || v === "0") continue; // missing
      const n = Number(v);
      if (Number.isFinite(n)) freqs[loci[l]!]![allele] = n;
    }
  }
  return { loci, freqs, N };
}

function compareCsvOutputs(
  dataset: string,
  pathPrefix: string,
  rText: string,
  tsText: string,
  /** When true, decode 3-char-encoded R allele keys back to TS form (haploid point-allele case). */
  decodeRAlleles: boolean = false,
): void {
  const r = parseCsv(rText);
  const ts = parseCsv(tsText);
  if (r.loci.length !== ts.loci.length) {
    recordMismatch({
      dataset,
      path: `${pathPrefix}.loci.length`,
      r: r.loci.length,
      ts: ts.loci.length,
      reason: "differs",
    });
  }
  for (const loc of new Set([...r.loci, ...ts.loci])) {
    if (r.N[loc] !== null && ts.N[loc] !== null) {
      compareNumbers(dataset, `${pathPrefix}.${loc}.N`, r.N[loc], ts.N[loc], DEFAULT_TOL);
    }
    const rFreq = r.freqs[loc] ?? {};
    const tsFreq = ts.freqs[loc] ?? {};
    if (decodeRAlleles) {
      // R-side keys are encoded ("90", "100", "930"); decode back to TS form
      // to align — but the encoding is lossy ("100" could be "10" or "10.0").
      // We additionally fall back to a sorted-multiset value comparison.
      const remapped: Record<string, number> = {};
      for (const k of Object.keys(rFreq)) {
        remapped[decodeHaploidAllele(k)] = rFreq[k]!;
      }
      const alleles = new Set([...Object.keys(remapped), ...Object.keys(tsFreq)]);
      for (const a of alleles) {
        // If only one side has it, also check whether the value matches anywhere
        // — accept silently if the multisets agree (handled below).
        if (remapped[a] !== undefined && tsFreq[a] !== undefined) {
          compareNumbers(dataset, `${pathPrefix}.${loc}.${a}`, remapped[a], tsFreq[a], DEFAULT_TOL);
        }
      }
      // Multiset value check: sorted frequencies should match.
      const rVals = Object.values(rFreq).sort((a, b) => a - b);
      const tsVals = Object.values(tsFreq).sort((a, b) => a - b);
      compareNumberArray(dataset, `${pathPrefix}.${loc}.values_sorted`, rVals, tsVals, DEFAULT_TOL);
    } else {
      const alleles = new Set([...Object.keys(rFreq), ...Object.keys(tsFreq)]);
      for (const a of alleles) {
        compareNumbers(dataset, `${pathPrefix}.${loc}.${a}`, rFreq[a], tsFreq[a], DEFAULT_TOL);
      }
    }
  }
}

function compareConversions(
  dataset: string,
  rConv: Record<string, unknown> | undefined,
  tsConv: Record<string, unknown> | undefined,
): void {
  if (!rConv || !tsConv) return;

  if (typeof rConv.genepop === "string" && typeof tsConv.genepop === "string") {
    compareGenepopOutputs(dataset, rConv.genepop, tsConv.genepop);
  }

  // Familias is per-population. Skip on haploid (R's straf2familias has a
  // diploid-only column-pairing assumption).
  if (rConv.familias && tsConv.familias) {
    const rFam = rConv.familias as Record<string, string | null>;
    const tsFam = tsConv.familias as Record<string, string | null>;
    for (const pop of new Set([...Object.keys(rFam), ...Object.keys(tsFam)])) {
      const rt = rFam[pop];
      const tt = tsFam[pop];
      if (typeof rt !== "string" || typeof tt !== "string") continue;
      compareFamiliasOutputs(dataset, `conversions.familias.${pop}`, rt, tt);
    }
  }

  for (const variant of ["euroformix", "strmix", "lrmix"] as const) {
    const rVar = rConv[variant] as Record<string, string | null> | undefined;
    const tsVar = tsConv[variant] as Record<string, string | null> | undefined;
    if (!rVar || !tsVar) continue;
    for (const pop of new Set([...Object.keys(rVar), ...Object.keys(tsVar)])) {
      const rt = rVar[pop];
      const tt = tsVar[pop];
      if (typeof rt !== "string" || typeof tt !== "string") continue;
      compareCsvOutputs(
        dataset,
        `conversions.${variant}.${pop}`,
        rt,
        tt,
        // For haploid point-allele datasets, R's allele keys went through
        // adegenet's lossy 3-char encoding; decode for comparison.
        isHaploidPointDataset(dataset),
      );
    }
  }

  // Arlequin: R's straf2arlequin has a known encoding bug — when no allele
  // in the locus contains a dot, it leaves alleles as raw strings ("10")
  // instead of 3-padded ("100"); when *any* allele contains a dot, it pads
  // to 2-3 chars inconsistently. TS produces consistent 3-char encoding
  // matching the Arlequin spec. We don't compare Arlequin byte-for-byte;
  // see the EXPECTED_DIFFERENCES note.
  // Mark every Arlequin path as expected-different so it's tracked.
  if (typeof rConv.arlequin === "string" && typeof tsConv.arlequin === "string") {
    if (rConv.arlequin !== tsConv.arlequin) {
      recordMismatch({
        dataset,
        path: "conversions.arlequin.text",
        r: `(R, ${(rConv.arlequin as string).length} chars)`,
        ts: `(TS, ${(tsConv.arlequin as string).length} chars)`,
        reason: "Arlequin output not byte-identical (R has variable-length allele encoding bug; TS uses spec-compliant 3-digit form).",
      });
    }
  }
}

// --------- Main ------------------------------------------------------------

const rFiles = (() => {
  try {
    return readdirSync(rDir).filter((f) => f.endsWith(".json"));
  } catch {
    return [];
  }
})();
const tsFiles = (() => {
  try {
    return readdirSync(tsDir).filter((f) => f.endsWith(".json"));
  } catch {
    return [];
  }
})();

if (rFiles.length === 0) {
  console.error(`No R-side dumps in ${rDir}.`);
  console.error(`Run \`Rscript straf-web/parity/dump_r.R\` from the repo root, then re-run this script.`);
  process.exit(2);
}
if (tsFiles.length === 0) {
  console.error(`No TS-side dumps in ${tsDir}.`);
  console.error(`Run \`node --import tsx parity/dump_ts.ts\` from straf-web/, then re-run this script.`);
  process.exit(2);
}

const all = new Set<string>([...rFiles, ...tsFiles]);
let datasetCount = 0;
const startTime = Date.now();
for (const f of all) {
  const dataset = f.replace(/\.json$/, "");
  if (!rFiles.includes(f)) {
    recordMismatch({ dataset, path: "(file)", r: "missing", ts: "present", reason: "no R dump" });
    continue;
  }
  if (!tsFiles.includes(f)) {
    recordMismatch({ dataset, path: "(file)", r: "present", ts: "missing", reason: "no TS dump" });
    continue;
  }
  datasetCount++;
  const r = JSON.parse(readFileSync(join(rDir, f), "utf8"));
  const ts = JSON.parse(readFileSync(join(tsDir, f), "utf8"));

  comparePerPopIndices(dataset, r.indices ?? {}, ts.indices ?? {});
  comparePerPopFreqs(dataset, r.freqs ?? {}, ts.freqs ?? {});
  if (r.perLocusFstats || ts.perLocusFstats) {
    comparePerLocusFstats(dataset, r.perLocusFstats ?? {}, ts.perLocusFstats ?? {});
  }
  if (r.pairwiseFst || ts.pairwiseFst) {
    comparePairwiseFst(dataset, r.pairwiseFst, ts.pairwiseFst);
  }
  if (r.distances || ts.distances) {
    compareDistances(dataset, r.distances ?? {}, ts.distances ?? {});
  }
  if (r.pca || ts.pca) {
    comparePca(dataset, r.pca, ts.pca);
  }
  if (r.conversions || ts.conversions) {
    compareConversions(dataset, r.conversions, ts.conversions);
  }
  if (r.haplotype || ts.haplotype) {
    compareHaplotype(dataset, r.haplotype ?? {}, ts.haplotype ?? {});
  }
}

classifyMismatches();
const unexpected = mismatches.filter((m) => !m.expected);
const expected = mismatches.filter((m) => m.expected);

const tookMs = Date.now() - startTime;
const showExpected = process.argv.includes("--show-expected");
console.log(`Compared ${datasetCount} dataset(s) in ${tookMs}ms.`);
console.log(
  `  unexpected mismatches: ${unexpected.length}` +
    (expected.length ? `   (plus ${expected.length} expected — see EXPECTED_DIFFERENCES)` : ""),
);

function printGroup(label: string, items: Mismatch[]): void {
  if (items.length === 0) return;
  console.log(`\n${label}:`);
  const byDataset = new Map<string, Mismatch[]>();
  for (const m of items) {
    if (!byDataset.has(m.dataset)) byDataset.set(m.dataset, []);
    byDataset.get(m.dataset)!.push(m);
  }
  for (const [ds, arr] of byDataset) {
    console.log(`\n  ${ds} (${arr.length}):`);
    const limit = 25;
    for (const m of arr.slice(0, limit)) {
      const detail =
        typeof m.r === "number" && typeof m.ts === "number"
          ? `r=${m.r}  ts=${m.ts}  diff=${m.diff?.toExponential(3)}`
          : `r=${JSON.stringify(m.r)}  ts=${JSON.stringify(m.ts)}`;
      const tag = m.expected ? `EXPECTED — ${m.expected}` : m.reason;
      console.log(`    ${m.path}: ${detail}  [${tag}]`);
    }
    if (arr.length > limit) console.log(`    ... and ${arr.length - limit} more`);
  }
}

printGroup("UNEXPECTED MISMATCHES", unexpected);
if (showExpected) {
  printGroup("EXPECTED MISMATCHES (filtered)", expected);
} else if (expected.length > 0) {
  // Just summarise per dataset + reason.
  const byReason = new Map<string, Map<string, number>>();
  for (const m of expected) {
    const key = m.expected!;
    if (!byReason.has(key)) byReason.set(key, new Map());
    const ds = byReason.get(key)!;
    ds.set(m.dataset, (ds.get(m.dataset) ?? 0) + 1);
  }
  console.log(`\nExpected mismatches (run with --show-expected for full detail):`);
  for (const [reason, ds] of byReason) {
    const total = Array.from(ds.values()).reduce((a, b) => a + b, 0);
    console.log(`  · [${total}] ${reason}`);
    for (const [d, n] of ds) console.log(`      ${d}: ${n}`);
  }
}

if (unexpected.length === 0) {
  console.log("\nAll checks passed (modulo expected differences).");
  process.exit(0);
}
process.exit(1);
