/**
 * Haplotype-level statistics for haploid (e.g. Y-STR / mtDNA) data.
 *
 * Mirrors `helpers_haplo_stats.R` from the original STRAF:
 *
 *   - haplotype frequency table (count + frequency, with H001-style labels)
 *   - haplotype diversity h = (n / (n−1)) · (1 − Σ p_i²)
 *     equivalently expressed in the original code as
 *     h = (n / (n−1)) · Σ p_i², which matches Nei's gene diversity for
 *     haplotypes (note the apparent sign — see source). We use the
 *     diversity convention 1 − Σ p_i² with the unbiased correction since
 *     that is what Nei (1987) recommends and what the report intent of
 *     the original tab implies. The raw Σ p_i² is reported alongside as
 *     "match probability" for traceability.
 *   - pairwise number of locus differences between distinct haplotypes
 *     (and proportion of differences = count / nLoci)
 */

import type { Genotypes } from "../parser.ts";
import { popMask, allMask } from "../parser.ts";

export interface HaplotypeRow {
  id: string;
  /** Underscore-joined alleles, with `null` rendered as "0". */
  haplotype: string;
  count: number;
  frequency: number;
}

export interface HaplotypeStats {
  /** Total number of complete haplotypes in the population. */
  n: number;
  nLoci: number;
  /** Distinct haplotypes (sorted by descending count). */
  haplotypes: HaplotypeRow[];
  /** Haplotype diversity h = n/(n−1) · (1 − Σ p²). */
  diversity: number;
  /** Σ p² — sum of squared haplotype frequencies (match probability). */
  matchProbability: number;
  /** pairDiff[i][j] = number of locus differences between haplotypes[i] and haplotypes[j]. */
  pairDiff: number[][];
  /** pairDiffProp[i][j] = pairDiff[i][j] / nLoci. */
  pairDiffProp: number[][];
}

export function haplotypeStats(genos: Genotypes, indMask?: boolean[]): HaplotypeStats {
  if (genos.ploidy !== 1) {
    throw new Error("Haplotype statistics require haploid data.");
  }

  const nLoci = genos.loci.length;
  const counts = new Map<string, number>();
  // Keep allele arrays for each unique haplotype so we can compute differences later.
  const haploAlleles = new Map<string, string[]>();
  let n = 0;

  for (let i = 0; i < genos.individuals.length; i++) {
    if (indMask && !indMask[i]) continue;
    const alleles: string[] = [];
    let complete = true;
    for (let l = 0; l < nLoci; l++) {
      const a = genos.alleles[l]![i]![0]!;
      if (a === null) {
        complete = false;
        break;
      }
      alleles.push(a);
    }
    if (!complete) continue;
    const key = alleles.join("|");
    counts.set(key, (counts.get(key) ?? 0) + 1);
    if (!haploAlleles.has(key)) haploAlleles.set(key, alleles);
    n++;
  }

  const sorted = Array.from(counts.entries()).sort((a, b) => b[1] - a[1]);
  const idWidth = String(sorted.length).length;
  const haplotypes: HaplotypeRow[] = sorted.map(([key, c], i) => ({
    id: `H${String(i + 1).padStart(idWidth, "0")}`,
    haplotype: key.split("|").join("|"),
    count: c,
    frequency: n > 0 ? c / n : 0,
  }));

  let sumP2 = 0;
  for (const h of haplotypes) sumP2 += h.frequency * h.frequency;
  const diversity = n > 1 ? (n / (n - 1)) * (1 - sumP2) : NaN;

  // Pairwise differences between distinct haplotypes.
  const k = haplotypes.length;
  const pairDiff: number[][] = Array.from({ length: k }, () => new Array(k).fill(0));
  const pairDiffProp: number[][] = Array.from({ length: k }, () => new Array(k).fill(0));
  for (let i = 0; i < k; i++) {
    const a = haploAlleles.get(sorted[i]![0])!;
    for (let j = i + 1; j < k; j++) {
      const b = haploAlleles.get(sorted[j]![0])!;
      let diff = 0;
      for (let l = 0; l < nLoci; l++) if (a[l] !== b[l]) diff++;
      pairDiff[i]![j] = diff;
      pairDiff[j]![i] = diff;
      pairDiffProp[i]![j] = diff / nLoci;
      pairDiffProp[j]![i] = diff / nLoci;
    }
  }

  return {
    n,
    nLoci,
    haplotypes,
    diversity,
    matchProbability: sumP2,
    pairDiff,
    pairDiffProp,
  };
}

export function haplotypeStatsForPop(
  genos: Genotypes,
  pop: string,
): HaplotypeStats {
  const mask = pop === "all" ? allMask(genos) : popMask(genos, pop);
  return haplotypeStats(genos, mask);
}
