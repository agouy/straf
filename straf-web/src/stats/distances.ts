import type { Genotypes } from "../parser.ts";
import { perLocusFrequencies } from "./freq.ts";

/**
 * Pairwise genetic distances between populations, computed from per-locus
 * allele frequencies. Returned distance matrix is square, symmetric, with
 * zeros on the diagonal.
 *
 * Methods (matching `adegenet::dist.genpop`):
 *
 *   - "nei"      Nei (1972) standard genetic distance
 *                D = -ln( Σ p_x p_y / sqrt(Σp_x² · Σp_y²) ), summed over loci
 *   - "rogers"   Rogers (1972) Euclidean-flavored distance
 *                D = (1/L) · Σ_l sqrt( ½ Σ_a (p_xla − p_yla)² )
 *   - "provesti" Provesti (1975) absolute genetic distance
 *                D = (1 / 2L) · Σ_l Σ_a |p_xla − p_yla|
 *   - "edwards"  Edwards (1971) angular distance
 *                D = sqrt( 1 − (1/L) · Σ_l Σ_a sqrt(p_xla · p_yla) )
 */

export type DistanceMethod = "nei" | "rogers" | "provesti" | "edwards" | "reynolds";

export interface PopulationDistances {
  populations: string[];
  /** matrix[i][j] = distance between populations[i] and populations[j] */
  matrix: number[][];
  method: DistanceMethod;
}

export function populationDistances(
  genos: Genotypes,
  method: DistanceMethod,
): PopulationDistances {
  const pops = Array.from(new Set(genos.populations));
  if (pops.length < 2) {
    throw new Error("At least 2 populations are required to compute distances.");
  }

  // Per-population, per-locus frequency maps.
  const freqsByPop: Map<string, number>[][] = pops.map((p) => {
    const mask = genos.populations.map((q) => q === p);
    return perLocusFrequencies(genos, mask).freqs;
  });

  const L = genos.loci.length;
  const matrix: number[][] = pops.map(() => new Array(pops.length).fill(0));

  for (let i = 0; i < pops.length; i++) {
    for (let j = i + 1; j < pops.length; j++) {
      const d = pairwiseDistance(freqsByPop[i]!, freqsByPop[j]!, L, method);
      matrix[i]![j] = d;
      matrix[j]![i] = d;
    }
  }

  return { populations: pops, matrix, method };
}

function pairwiseDistance(
  Fx: Map<string, number>[],
  Fy: Map<string, number>[],
  L: number,
  method: DistanceMethod,
): number {
  switch (method) {
    case "nei": {
      // Σ over loci of p_x·p_y, p_x², p_y²
      let num = 0;
      let xx = 0;
      let yy = 0;
      let usedLoci = 0;
      for (let l = 0; l < L; l++) {
        const fx = Fx[l]!;
        const fy = Fy[l]!;
        if (fx.size === 0 || fy.size === 0) continue;
        const alleles = new Set([...fx.keys(), ...fy.keys()]);
        for (const a of alleles) {
          const px = fx.get(a) ?? 0;
          const py = fy.get(a) ?? 0;
          num += px * py;
          xx += px * px;
          yy += py * py;
        }
        usedLoci++;
      }
      if (usedLoci === 0 || xx === 0 || yy === 0 || num === 0) return Infinity;
      const I = num / Math.sqrt(xx * yy);
      return -Math.log(Math.min(1, Math.max(1e-12, I)));
    }
    case "rogers": {
      let sumLoc = 0;
      let usedLoci = 0;
      for (let l = 0; l < L; l++) {
        const fx = Fx[l]!;
        const fy = Fy[l]!;
        if (fx.size === 0 || fy.size === 0) continue;
        const alleles = new Set([...fx.keys(), ...fy.keys()]);
        let s = 0;
        for (const a of alleles) {
          const px = fx.get(a) ?? 0;
          const py = fy.get(a) ?? 0;
          s += (px - py) * (px - py);
        }
        sumLoc += Math.sqrt(0.5 * s);
        usedLoci++;
      }
      return usedLoci === 0 ? 0 : sumLoc / usedLoci;
    }
    case "provesti": {
      let sumLoc = 0;
      let usedLoci = 0;
      for (let l = 0; l < L; l++) {
        const fx = Fx[l]!;
        const fy = Fy[l]!;
        if (fx.size === 0 || fy.size === 0) continue;
        const alleles = new Set([...fx.keys(), ...fy.keys()]);
        for (const a of alleles) {
          const px = fx.get(a) ?? 0;
          const py = fy.get(a) ?? 0;
          sumLoc += Math.abs(px - py);
        }
        usedLoci++;
      }
      return usedLoci === 0 ? 0 : sumLoc / (2 * usedLoci);
    }
    case "edwards": {
      let sumLoc = 0;
      let usedLoci = 0;
      for (let l = 0; l < L; l++) {
        const fx = Fx[l]!;
        const fy = Fy[l]!;
        if (fx.size === 0 || fy.size === 0) continue;
        const alleles = new Set([...fx.keys(), ...fy.keys()]);
        let s = 0;
        for (const a of alleles) {
          const px = fx.get(a) ?? 0;
          const py = fy.get(a) ?? 0;
          s += Math.sqrt(px * py);
        }
        sumLoc += s;
        usedLoci++;
      }
      if (usedLoci === 0) return 0;
      return Math.sqrt(Math.max(0, 1 - sumLoc / usedLoci));
    }
    case "reynolds": {
      // Reynolds, Weir & Cockerham (1983) coancestrality coefficient θ_W.
      // Per locus: num = ½ Σ_a (p_xa − p_ya)²
      //           den = 1 − Σ_a p_xa · p_ya
      // Distance = − ln( 1 − Σ num / Σ den ) summed over loci.
      let num = 0;
      let den = 0;
      let usedLoci = 0;
      for (let l = 0; l < L; l++) {
        const fx = Fx[l]!;
        const fy = Fy[l]!;
        if (fx.size === 0 || fy.size === 0) continue;
        const alleles = new Set([...fx.keys(), ...fy.keys()]);
        let nLoc = 0;
        let dLoc = 1;
        for (const a of alleles) {
          const px = fx.get(a) ?? 0;
          const py = fy.get(a) ?? 0;
          nLoc += (px - py) * (px - py);
          dLoc -= px * py;
        }
        num += nLoc / 2;
        den += dLoc;
        usedLoci++;
      }
      if (usedLoci === 0 || den <= 0) return 0;
      const ratio = num / den;
      if (ratio >= 1) return Infinity;
      return -Math.log(Math.max(1e-12, 1 - ratio));
    }
  }
}
