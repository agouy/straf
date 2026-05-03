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
