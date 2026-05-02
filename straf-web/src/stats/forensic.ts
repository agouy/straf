import type { Genotypes } from "../parser.ts";
import { perLocusFrequencies } from "./freq.ts";

/**
 * Forensic and population genetics indices, per locus.
 *
 * Formulas (matching pegas / adegenet / hierfstat as used by STRAF):
 *
 *   N        = number of typed gene copies at the locus (sum over individuals
 *              of non-missing allele counts)
 *   Nall     = number of distinct alleles observed
 *   GD (Hexp)= 1 - sum(p_i^2), with sample-size correction N/(N-1)
 *              (Nei's unbiased gene diversity / expected heterozygosity)
 *   PIC      = polymorphism information content
 *              = 1 - sum(p_i^2) - sum_{i<j} 2 * p_i^2 * p_j^2
 *   PM       = matching probability = sum_g (g_count/N_ind)^2 over genotypes
 *   PD       = power of discrimination = 1 - PM
 *   Hobs     = observed heterozygosity (diploid only)
 *   PE       = power of exclusion (diploid only)
 *              = Hobs^2 * (1 - 2 * Hobs * (1 - Hobs)^2)
 *   TPI      = typical paternity index (diploid only) = 1 / (2 * (1 - Hobs))
 */
export interface LocusStats {
  locus: string;
  N: number;
  Nall: number;
  GD: number; // a.k.a. Hexp (unbiased)
  PIC: number;
  PM: number;
  PD: number;
  Hobs?: number;
  PE?: number;
  TPI?: number;
}

export function locusIndices(genos: Genotypes, indMask?: boolean[]): LocusStats[] {
  const { freqs, N } = perLocusFrequencies(genos, indMask);
  const out: LocusStats[] = [];

  for (let l = 0; l < genos.loci.length; l++) {
    const fmap = freqs[l]!;
    const Nl = N[l]!;
    const ps = Array.from(fmap.values());
    const sumP2 = ps.reduce((s, p) => s + p * p, 0);

    // PIC
    let sum2pq2 = 0;
    for (let i = 0; i < ps.length; i++) {
      for (let j = i + 1; j < ps.length; j++) {
        sum2pq2 += 2 * ps[i]! * ps[i]! * ps[j]! * ps[j]!;
      }
    }
    const PIC = 1 - sumP2 - sum2pq2;

    // GD with sample-size correction
    const GDraw = 1 - sumP2;
    const GD = Nl > 1 ? (GDraw * Nl) / (Nl - 1) : GDraw;

    // PM and PD: enumerate observed unordered genotypes
    const tuples = genos.alleles[l]!;
    const genoCounts = new Map<string, number>();
    let nGeno = 0;
    let hetCount = 0;
    let nDiploidTyped = 0;
    for (let i = 0; i < tuples.length; i++) {
      if (indMask && !indMask[i]) continue;
      const t = tuples[i]!;
      if (t.some((a) => a === null)) continue; // skip incomplete genotypes
      const sorted = [...(t as string[])].sort();
      const key = sorted.join("/");
      genoCounts.set(key, (genoCounts.get(key) ?? 0) + 1);
      nGeno++;
      if (genos.ploidy === 2) {
        nDiploidTyped++;
        if (sorted[0] !== sorted[1]) hetCount++;
      }
    }
    let PM = 0;
    if (nGeno > 0) {
      for (const c of genoCounts.values()) PM += (c / nGeno) ** 2;
    }
    const PD = 1 - PM;

    const stat: LocusStats = {
      locus: genos.loci[l]!,
      N: Nl,
      Nall: ps.length,
      GD,
      PIC,
      PM,
      PD,
    };

    if (genos.ploidy === 2 && nDiploidTyped > 0) {
      const Hobs = hetCount / nDiploidTyped;
      stat.Hobs = Hobs;
      stat.PE = Hobs * Hobs * (1 - 2 * Hobs * (1 - Hobs) * (1 - Hobs));
      stat.TPI = Hobs < 1 ? 1 / (2 * (1 - Hobs)) : Infinity;
    }

    out.push(stat);
  }

  return out;
}

/** Combined forensic indices over all loci (cumulative PD and PE). */
export function combinedIndices(stats: LocusStats[]): {
  combinedPM: number;
  combinedPD: number;
  combinedPE?: number;
} {
  let prodPM = 1;
  for (const s of stats) prodPM *= s.PM;
  const combinedPD = 1 - prodPM;

  let combinedPE: number | undefined;
  if (stats.every((s) => s.PE !== undefined)) {
    let prod = 1;
    for (const s of stats) prod *= 1 - (s.PE as number);
    combinedPE = 1 - prod;
  }
  return { combinedPM: prodPM, combinedPD, combinedPE };
}
