/**
 * Parity dump (TS side).
 *
 * Mirrors `dump_r.R` — runs the new TypeScript implementations on every
 * dataset listed in parity/manifest.json and writes one JSON file per dataset
 * under parity/results_ts/, in a layout the comparator can match against the
 * R-side dumps.
 *
 * Run from straf-web/:
 *   node --import tsx parity/dump_ts.ts
 */

import { readFileSync, writeFileSync, mkdirSync, readdirSync } from "node:fs";
import { join, dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

import { parseStrafTsv, allMask, popMask, uniquePopulations } from "../src/parser.ts";
import type { Genotypes, Ploidy } from "../src/parser.ts";
import { alleleFrequencies } from "../src/stats/freq.ts";
import { locusIndices } from "../src/stats/forensic.ts";
import { perLocusFStats, pairwiseFst, basicStatsPerLocus } from "../src/stats/fst.ts";
import { populationDistances } from "../src/stats/distances.ts";
import type { DistanceMethod } from "../src/stats/distances.ts";
import { pca } from "../src/stats/pca.ts";
import { haplotypeStatsForPop } from "../src/stats/haplotype.ts";
import {
  toGenepop,
  toArlequin,
  toFamilias,
  toAlleleFreqCsv,
} from "../src/convert/formats.ts";

const here = dirname(fileURLToPath(import.meta.url));
const dsDir = join(here, "datasets");
const outDir = join(here, "results_ts");
mkdirSync(outDir, { recursive: true });

interface ManifestEntry {
  name: string;
  ploidy: 1 | 2;
  description: string;
}
const manifest: { datasets: ManifestEntry[] } = JSON.parse(
  readFileSync(join(here, "manifest.json"), "utf8"),
);

function maskFor(g: Genotypes, pop: string): boolean[] {
  return pop === "all" ? allMask(g) : popMask(g, pop);
}

/** Per-population indices, keyed { pop: { locus: { N, Nall, GD, PIC, ... } } }. */
function dumpIndices(g: Genotypes): Record<string, Record<string, Record<string, number | null>>> {
  const out: Record<string, Record<string, Record<string, number | null>>> = {};
  const pops = ["all", ...uniquePopulations(g)];
  for (const p of pops) {
    const stats = locusIndices(g, maskFor(g, p));
    const byLocus: Record<string, Record<string, number | null>> = {};
    for (const s of stats) {
      const row: Record<string, number | null> = {
        N: s.N,
        Nall: s.Nall,
        GD: s.GD,
        PIC: s.PIC,
        PM: s.PM,
        PD: s.PD,
      };
      if (s.Hobs !== undefined) row.Hobs = s.Hobs;
      if (s.PE !== undefined) row.PE = s.PE;
      if (s.TPI !== undefined) row.TPI = Number.isFinite(s.TPI) ? s.TPI : null;
      byLocus[s.locus] = row;
    }
    out[p] = byLocus;
  }
  return out;
}

/** Per-population frequency tables, keyed { pop: { locus: { N, freq: { allele: p } } } }. */
function dumpFreqs(g: Genotypes): Record<string, Record<string, { N: number; freq: Record<string, number> }>> {
  const out: Record<string, Record<string, { N: number; freq: Record<string, number> }>> = {};
  const pops = ["all", ...uniquePopulations(g)];
  for (const p of pops) {
    const tbl = alleleFrequencies(g, maskFor(g, p));
    const byLocus: Record<string, { N: number; freq: Record<string, number> }> = {};
    for (let l = 0; l < tbl.loci.length; l++) {
      const freq: Record<string, number> = {};
      for (let a = 0; a < tbl.alleles.length; a++) {
        const v = tbl.values[a]![l];
        if (v !== null && v !== undefined) freq[tbl.alleles[a]!] = v;
      }
      byLocus[tbl.loci[l]!] = { N: tbl.N[l]!, freq };
    }
    out[p] = byLocus;
  }
  return out;
}

function dumpPerLocusFstats(g: Genotypes): Record<string, Record<string, number | null>> {
  // Mirror the R STRAF UI: Ht/Hs/Fis come from basic.stats; Fst from W&C wc().
  const wc = perLocusFStats(g);
  const bs = basicStatsPerLocus(g);
  const out: Record<string, Record<string, number | null>> = {};
  const wcByLoc = new Map(wc.map((f) => [f.locus, f] as const));
  for (const b of bs) {
    const w = wcByLoc.get(b.locus);
    out[b.locus] = {
      Ho: numOrNull(b.Ho),
      Hs: numOrNull(b.Hs),
      Ht: numOrNull(b.Ht),
      Fis: numOrNull(b.Fis),
      Fst: w ? numOrNull(w.Fst) : null,
    };
  }
  return out;
}

function dumpPairwiseFst(g: Genotypes): { populations: string[]; matrix: (number | null)[][] } {
  const r = pairwiseFst(g);
  return {
    populations: r.populations,
    matrix: r.matrix.map((row) => row.map(numOrNull)),
  };
}

function dumpDist(g: Genotypes, method: DistanceMethod): { populations: string[]; matrix: number[][] } {
  const d = populationDistances(g, method);
  return { populations: d.populations, matrix: d.matrix };
}

/**
 * Dump PCA in a sign/rotation-invariant form: variance proportions,
 * eigenvalues, and a pairwise-distance matrix on the first up-to-5 PCs.
 *
 * The original (R) PCA runs `prcomp(scale=TRUE)` on `adegenet::makefreq(...)`
 * (frequencies in [0, 1] with mean imputation). Our TS PCA builds a *dosage*
 * matrix (counts 0/1/2) but uses the same mean imputation and centering, so:
 *   - variance proportions: identical
 *   - eigenvalues: differ by a constant factor when scale=FALSE; when scale=TRUE
 *     (the original's choice), the per-column standardisation removes that
 *     factor and eigenvalues should match too
 *   - pairwise distances on PC space: match exactly under scale=TRUE
 */
function dumpPca(g: Genotypes): {
  nInd: number;
  nCols: number;
  eigenvalues: number[];
  variance: number[];
  cumulative: number[];
  pairDist: number[][];
  individuals: string[];
} {
  const res = pca(g, undefined, { scale: true });
  // Adegenet's PCA passes a scaled frequency matrix to prcomp, so eigenvalues
  // are on the same scale as λ_i = sdev_i^2 = (variance per column-scaled PC).
  // Our `scale=true` variant matches that.
  // Eigenvalues from our Jacobi are sums-of-squares (i.e. (n-1) * variance).
  // Convert to "variance" by dividing by (n-1) so they line up with prcomp.
  const denom = Math.max(1, res.individuals.length - 1);
  const eigenvalues = res.eigenvalues.map((v) => v / denom);
  const k = Math.min(5, res.scores[0]!.length);
  const n = res.scores.length;
  const pairDist: number[][] = Array.from({ length: n }, () => new Array(n).fill(0) as number[]);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      let s = 0;
      for (let c = 0; c < k; c++) {
        const d = res.scores[i]![c]! - res.scores[j]![c]!;
        s += d * d;
      }
      const dij = Math.sqrt(s);
      pairDist[i]![j] = dij;
      pairDist[j]![i] = dij;
    }
  }
  let acc = 0;
  const cumulative = res.variance.map((v) => (acc += v));
  return {
    nInd: res.individuals.length,
    nCols: res.alleleLabels.length,
    eigenvalues,
    variance: res.variance,
    cumulative,
    pairDist,
    individuals: res.individuals,
  };
}

function dumpHaplotype(
  g: Genotypes,
): Record<
  string,
  {
    n: number;
    nLoci: number;
    haplotypes: { id: string; haplotype: string; count: number; frequency: number }[];
    sumP2: number;
    h_div_R: number;
    diversity_TSlike: number;
    d_mat: number[][];
  }
> {
  const out: Record<string, ReturnType<typeof packHaploEntry>> = {};
  const pops = ["all", ...uniquePopulations(g)];
  for (const p of pops) {
    const r = haplotypeStatsForPop(g, p);
    out[p] = packHaploEntry(r.n, r.nLoci, r.haplotypes, r.matchProbability, r.diversity, r.pairDiff);
  }
  return out;
}

function packHaploEntry(
  n: number,
  nLoci: number,
  haps: { id: string; haplotype: string; count: number; frequency: number }[],
  sumP2: number,
  diversityTsLike: number,
  d_mat: number[][],
) {
  return {
    n,
    nLoci,
    haplotypes: haps.map((h) => ({ id: h.id, haplotype: h.haplotype, count: h.count, frequency: h.frequency })),
    sumP2,
    // The R version computes h_div = n/(n-1) * sum(p^2) (i.e. a sample-corrected
    // *match probability*, not Nei's gene diversity). Reproduce that here so
    // we can compare apples to apples.
    h_div_R: n > 1 ? (n / (n - 1)) * sumP2 : 0,
    // Our TS convention: Nei's unbiased gene diversity n/(n-1) * (1 - Σp²).
    diversity_TSlike: diversityTsLike,
    d_mat,
  };
}

function numOrNull(v: number): number | null {
  return Number.isFinite(v) ? v : null;
}

// --------- Run --------------------------------------------------------------

const datasets = manifest.datasets;
let count = 0;
for (const ds of datasets) {
  const path = join(dsDir, `${ds.name}.txt`);
  const text = readFileSync(path, "utf8");
  const ploidy = ds.ploidy as Ploidy;
  process.stderr.write(`[TS] ${ds.name} (ploidy=${ploidy})\n`);
  const g = parseStrafTsv(text, ploidy);

  const out: Record<string, unknown> = {
    dataset: ds.name,
    ploidy,
    nInd: g.individuals.length,
    pops: uniquePopulations(g).sort(),
    indices: dumpIndices(g),
    freqs: dumpFreqs(g),
  };

  const npop = uniquePopulations(g).length;
  if (ploidy === 2 && npop > 1 && g.loci.length > 1) {
    try {
      out.perLocusFstats = dumpPerLocusFstats(g);
    } catch {
      out.perLocusFstats = null;
    }
    try {
      out.pairwiseFst = dumpPairwiseFst(g);
    } catch {
      out.pairwiseFst = null;
    }
    out.distances = {};
    for (const m of ["nei", "edwards", "reynolds", "rogers", "provesti"] as const) {
      try {
        (out.distances as Record<string, unknown>)[m] = dumpDist(g, m);
      } catch {
        (out.distances as Record<string, unknown>)[m] = null;
      }
    }
  }

  if (g.individuals.length >= 2) {
    try {
      out.pca = dumpPca(g);
    } catch (err) {
      out.pca = null;
      process.stderr.write(`  PCA skipped: ${(err as Error).message}\n`);
    }
  }

  if (ploidy === 1) {
    try {
      out.haplotype = dumpHaplotype(g);
    } catch (err) {
      out.haplotype = null;
      process.stderr.write(`  haplotype skipped: ${(err as Error).message}\n`);
    }
  }

  // File conversions: capture the raw text of each output. Skips match
  // dump_r.R: Arlequin and Familias are diploid-only on the R side.
  const conversions: Record<string, unknown> = {};
  try {
    conversions.genepop = toGenepop(g);
  } catch {
    conversions.genepop = null;
  }
  if (ploidy === 2) {
    try {
      conversions.arlequin = toArlequin(g);
    } catch {
      conversions.arlequin = null;
    }
    const familias: Record<string, string | null> = {};
    for (const p of ["all", ...uniquePopulations(g)]) {
      try {
        familias[p] = toFamilias(g, p);
      } catch {
        familias[p] = null;
      }
    }
    conversions.familias = familias;
  }
  for (const variant of ["euroformix", "strmix", "lrmix"] as const) {
    const tbl: Record<string, string | null> = {};
    for (const p of ["all", ...uniquePopulations(g)]) {
      try {
        tbl[p] = toAlleleFreqCsv(g, variant, p);
      } catch {
        tbl[p] = null;
      }
    }
    conversions[variant] = tbl;
  }
  out.conversions = conversions;

  writeFileSync(join(outDir, `${ds.name}.json`), JSON.stringify(out, null, 2));
  count++;
}

process.stderr.write(`[TS] wrote ${count} JSON files to ${resolve(outDir)}\n`);

// Friendly: if R hasn't dumped yet, suggest it.
const rDir = join(here, "results_r");
const rCount = (() => {
  try {
    return readdirSync(rDir).filter((f) => f.endsWith(".json")).length;
  } catch {
    return 0;
  }
})();
if (rCount === 0) {
  process.stderr.write(`[TS] note: no R-side dumps found in ${rDir}. Run \`Rscript straf-web/parity/dump_r.R\` from the repo root.\n`);
}
