import type { Genotypes } from "../parser.ts";

/**
 * Pairwise Weir & Cockerham (1984) θ̂ between every pair of populations.
 *
 * For r populations with sample sizes n_1, ..., n_r, per locus and per allele:
 *
 *   n̄ = mean(n_i)
 *   n_c = (Σ n_i − Σ n_i² / Σ n_i) / (r − 1)
 *   p̃_i  = sample allele frequency in population i  (count / 2 n_i for diploid)
 *   p̄    = (Σ n_i p̃_i) / (Σ n_i)
 *   h̄    = (Σ n_i h_i) / (Σ n_i)   where h_i = fraction of individuals
 *                                  heterozygous for the allele in pop i
 *   s²   = (Σ n_i (p̃_i − p̄)²) / ((r − 1) n̄)
 *
 *   a = n̄/n_c · ( s² − 1/(n̄−1) · ( p̄(1−p̄) − (r−1)/r · s² − h̄/4 ) )
 *   b = n̄/(n̄−1) · ( p̄(1−p̄) − (r−1)/r · s² − (2 n̄ − 1)/(4 n̄) · h̄ )
 *   c = h̄ / 2
 *
 * Then θ̂ = (Σ_loci Σ_alleles a) / (Σ_loci Σ_alleles (a + b + c)).
 *
 * For pairwise Fst we set r = 2 and restrict the populations to the pair.
 * Negative values (which can occur with small samples) are reported as-is;
 * many implementations clamp to zero — we don't, to match hierfstat output.
 */

export interface PairwiseFst {
  populations: string[];
  /** matrix[i][j] = Fst between populations[i] and populations[j]; diagonal = 0 */
  matrix: number[][];
}

export function pairwiseFst(genos: Genotypes): PairwiseFst {
  if (genos.ploidy !== 2) {
    throw new Error("Pairwise Fst (Weir & Cockerham) requires diploid data.");
  }
  const pops = Array.from(new Set(genos.populations));
  if (pops.length < 2) {
    throw new Error("At least 2 populations are required to compute pairwise Fst.");
  }

  const matrix: number[][] = pops.map(() => new Array(pops.length).fill(0));

  for (let i = 0; i < pops.length; i++) {
    for (let j = i + 1; j < pops.length; j++) {
      const f = wcTwoPop(genos, pops[i]!, pops[j]!);
      matrix[i]![j] = f;
      matrix[j]![i] = f;
    }
  }

  return { populations: pops, matrix };
}

function wcTwoPop(genos: Genotypes, popA: string, popB: string): number {
  const r = 2;
  let sumA = 0; // numerator
  let sumABC = 0; // denominator

  for (let l = 0; l < genos.loci.length; l++) {
    // Compute n_i (number of individuals typed at locus l in each pop), and
    // for each allele its count and the count of heterozygotes.
    const tuples = genos.alleles[l]!;
    let nA = 0;
    let nB = 0;
    const alleleA = new Map<string, number>();
    const alleleB = new Map<string, number>();
    const hetA = new Map<string, number>();
    const hetB = new Map<string, number>();

    for (let i = 0; i < tuples.length; i++) {
      const pop = genos.populations[i]!;
      const isA = pop === popA;
      const isB = pop === popB;
      if (!isA && !isB) continue;
      const t = tuples[i]!;
      if (t[0] === null || t[1] === null) continue;
      const a0 = t[0] as string;
      const a1 = t[1] as string;
      const aMap = isA ? alleleA : alleleB;
      const hMap = isA ? hetA : hetB;
      if (isA) nA++;
      else nB++;
      aMap.set(a0, (aMap.get(a0) ?? 0) + 1);
      aMap.set(a1, (aMap.get(a1) ?? 0) + 1);
      if (a0 !== a1) {
        hMap.set(a0, (hMap.get(a0) ?? 0) + 1);
        hMap.set(a1, (hMap.get(a1) ?? 0) + 1);
      }
    }

    if (nA === 0 || nB === 0) continue;

    const totalN = nA + nB;
    const nBar = totalN / r;
    const nC = (totalN - (nA * nA + nB * nB) / totalN) / (r - 1);
    if (nBar <= 1 || nC <= 0) continue;

    const allAlleles = new Set<string>();
    for (const a of alleleA.keys()) allAlleles.add(a);
    for (const a of alleleB.keys()) allAlleles.add(a);

    for (const a of allAlleles) {
      const c0 = alleleA.get(a) ?? 0;
      const c1 = alleleB.get(a) ?? 0;
      const p0 = c0 / (2 * nA);
      const p1 = c1 / (2 * nB);
      const pBar = (nA * p0 + nB * p1) / totalN;
      const h0 = (hetA.get(a) ?? 0) / nA;
      const h1 = (hetB.get(a) ?? 0) / nB;
      const hBar = (nA * h0 + nB * h1) / totalN;
      const s2 =
        (nA * (p0 - pBar) * (p0 - pBar) + nB * (p1 - pBar) * (p1 - pBar)) /
        ((r - 1) * nBar);

      const inner = pBar * (1 - pBar) - ((r - 1) / r) * s2;
      const aT = (nBar / nC) * (s2 - (1 / (nBar - 1)) * (inner - hBar / 4));
      const bT = (nBar / (nBar - 1)) * (inner - ((2 * nBar - 1) / (4 * nBar)) * hBar);
      const cT = hBar / 2;

      sumA += aT;
      sumABC += aT + bT + cT;
    }
  }

  if (sumABC === 0) return 0;
  return sumA / sumABC;
}

/**
 * Per-locus Weir & Cockerham F-statistics across all r populations.
 *
 * Returns Ht (total heterozygosity), Hs (within-pop heterozygosity),
 * Fis, Fst (= θ̂), Fit per locus.
 *
 * Definitions follow Nei (1987) / Weir & Cockerham (1984):
 *   Ht  = 1 − Σ p̄_a²        with p̄_a the weighted mean allele frequency
 *   Hs  = mean over pops of (1 − Σ p_a²) — sample-size weighted
 *   Fst = (a) / (a + b + c)  per locus, sum over alleles
 *   Fis = (b) / (b + c)
 *   Fit = (a + b) / (a + b + c)
 */
export interface LocusFStats {
  locus: string;
  Ht: number;
  Hs: number;
  Fis: number;
  Fst: number;
  Fit: number;
}

export function perLocusFStats(genos: Genotypes): LocusFStats[] {
  if (genos.ploidy !== 2) {
    throw new Error("Per-locus F-statistics require diploid data.");
  }
  const pops = Array.from(new Set(genos.populations));
  if (pops.length < 2) {
    throw new Error("At least 2 populations are required to compute F-statistics.");
  }
  const r = pops.length;
  const out: LocusFStats[] = [];

  for (let l = 0; l < genos.loci.length; l++) {
    const tuples = genos.alleles[l]!;

    const nByPop = new Map<string, number>();
    const alleleByPop = new Map<string, Map<string, number>>();
    const hetByPop = new Map<string, Map<string, number>>();

    for (const p of pops) {
      nByPop.set(p, 0);
      alleleByPop.set(p, new Map());
      hetByPop.set(p, new Map());
    }

    for (let i = 0; i < tuples.length; i++) {
      const pop = genos.populations[i]!;
      const t = tuples[i]!;
      if (t[0] === null || t[1] === null) continue;
      nByPop.set(pop, nByPop.get(pop)! + 1);
      const aMap = alleleByPop.get(pop)!;
      const hMap = hetByPop.get(pop)!;
      const a0 = t[0] as string;
      const a1 = t[1] as string;
      aMap.set(a0, (aMap.get(a0) ?? 0) + 1);
      aMap.set(a1, (aMap.get(a1) ?? 0) + 1);
      if (a0 !== a1) {
        hMap.set(a0, (hMap.get(a0) ?? 0) + 1);
        hMap.set(a1, (hMap.get(a1) ?? 0) + 1);
      }
    }

    let totalN = 0;
    for (const p of pops) totalN += nByPop.get(p)!;
    if (totalN === 0) {
      out.push({ locus: genos.loci[l]!, Ht: NaN, Hs: NaN, Fis: NaN, Fst: NaN, Fit: NaN });
      continue;
    }
    const nBar = totalN / r;
    let sumNi2 = 0;
    for (const p of pops) {
      const ni = nByPop.get(p)!;
      sumNi2 += ni * ni;
    }
    const nC = (totalN - sumNi2 / totalN) / (r - 1);

    const allAlleles = new Set<string>();
    for (const m of alleleByPop.values()) for (const a of m.keys()) allAlleles.add(a);

    let sumA = 0;
    let sumB = 0;
    let sumC = 0;
    let pSqMean = 0; // Σ_a p̄_a² — for Ht
    let HsAccum = 0; // Σ_pop n_i (1 − Σ_a p_ia²) / Σ n_i — for Hs

    for (const p of pops) {
      const ni = nByPop.get(p)!;
      if (ni === 0) continue;
      const aMap = alleleByPop.get(p)!;
      let sumP2 = 0;
      for (const c of aMap.values()) sumP2 += (c / (2 * ni)) ** 2;
      HsAccum += (ni / totalN) * (1 - sumP2);
    }

    for (const a of allAlleles) {
      // Build per-pop frequencies
      const ps: number[] = [];
      const hs: number[] = [];
      const ns: number[] = [];
      for (const p of pops) {
        const ni = nByPop.get(p)!;
        ns.push(ni);
        if (ni === 0) {
          ps.push(0);
          hs.push(0);
          continue;
        }
        const ca = alleleByPop.get(p)!.get(a) ?? 0;
        const ha = hetByPop.get(p)!.get(a) ?? 0;
        ps.push(ca / (2 * ni));
        hs.push(ha / ni);
      }
      let pBar = 0;
      let hBar = 0;
      for (let k = 0; k < r; k++) {
        pBar += (ns[k]! * ps[k]!) / totalN;
        hBar += (ns[k]! * hs[k]!) / totalN;
      }
      let s2 = 0;
      for (let k = 0; k < r; k++) {
        s2 += (ns[k]! * (ps[k]! - pBar) * (ps[k]! - pBar)) / ((r - 1) * nBar);
      }

      pSqMean += pBar * pBar;

      if (nBar <= 1 || nC <= 0) continue;
      const inner = pBar * (1 - pBar) - ((r - 1) / r) * s2;
      const aT = (nBar / nC) * (s2 - (1 / (nBar - 1)) * (inner - hBar / 4));
      const bT = (nBar / (nBar - 1)) * (inner - ((2 * nBar - 1) / (4 * nBar)) * hBar);
      const cT = hBar / 2;
      sumA += aT;
      sumB += bT;
      sumC += cT;
    }

    const denom = sumA + sumB + sumC;
    const Fst = denom !== 0 ? sumA / denom : NaN;
    const Fis = sumB + sumC !== 0 ? sumB / (sumB + sumC) : NaN;
    const Fit = denom !== 0 ? (sumA + sumB) / denom : NaN;
    const Ht = 1 - pSqMean;
    const Hs = HsAccum;

    out.push({
      locus: genos.loci[l]!,
      Ht,
      Hs,
      Fis,
      Fst,
      Fit,
    });
  }

  return out;
}
