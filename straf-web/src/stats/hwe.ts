import type { Genotypes } from "../parser.ts";
import { perLocusFrequencies } from "./freq.ts";

/**
 * Asymptotic chi-square Hardy–Weinberg equilibrium test, per locus.
 *
 * Compares observed unordered diploid genotype counts against expected counts
 * under HWE (p_i^2 for homozygotes, 2 p_i p_j for heterozygotes), using
 * sample frequencies as estimates of allele frequencies.
 *
 * Degrees of freedom = K(K-1)/2 where K is the number of distinct alleles.
 *
 * Notes:
 *   - This is an asymptotic approximation. The original STRAF uses Genepop's
 *     exact MCMC test (Guo & Thompson 1992), which is more accurate when
 *     expected genotype counts are small. Reimplementing the MCMC sampler
 *     would be a couple hundred lines and is left as a follow-up.
 *   - Genotypes with missing alleles are dropped from the test for that locus.
 *   - Returned p-value is NaN when df <= 0 or no complete genotypes exist.
 */
export interface HweResult {
  locus: string;
  N: number; // number of complete diploid genotypes
  K: number; // number of distinct alleles
  df: number;
  chiSq: number;
  pValue: number;
}

export function hweChiSquare(genos: Genotypes, indMask?: boolean[]): HweResult[] {
  if (genos.ploidy !== 2) {
    throw new Error("HWE test requires diploid data.");
  }

  const { freqs } = perLocusFrequencies(genos, indMask);
  const out: HweResult[] = [];

  for (let l = 0; l < genos.loci.length; l++) {
    const fmap = freqs[l]!;
    const alleles = Array.from(fmap.keys());
    const K = alleles.length;
    const tuples = genos.alleles[l]!;

    // Observed unordered genotype counts among complete diploid genotypes.
    const obs = new Map<string, number>();
    let N = 0;
    for (let i = 0; i < tuples.length; i++) {
      if (indMask && !indMask[i]) continue;
      const t = tuples[i]!;
      if (t[0] === null || t[1] === null) continue;
      const key = sortedKey(t[0] as string, t[1] as string);
      obs.set(key, (obs.get(key) ?? 0) + 1);
      N++;
    }

    const df = (K * (K - 1)) / 2;
    if (N === 0 || df <= 0) {
      out.push({ locus: genos.loci[l]!, N, K, df, chiSq: NaN, pValue: NaN });
      continue;
    }

    // Expected counts under HWE.
    let chiSq = 0;
    for (let i = 0; i < K; i++) {
      const ai = alleles[i]!;
      const pi = fmap.get(ai)!;
      // homozygote
      const eHom = N * pi * pi;
      const oHom = obs.get(sortedKey(ai, ai)) ?? 0;
      if (eHom > 0) chiSq += ((oHom - eHom) ** 2) / eHom;
      // heterozygotes
      for (let j = i + 1; j < K; j++) {
        const aj = alleles[j]!;
        const pj = fmap.get(aj)!;
        const eHet = N * 2 * pi * pj;
        const oHet = obs.get(sortedKey(ai, aj)) ?? 0;
        if (eHet > 0) chiSq += ((oHet - eHet) ** 2) / eHet;
      }
    }

    out.push({
      locus: genos.loci[l]!,
      N,
      K,
      df,
      chiSq,
      pValue: chiSquarePValue(chiSq, df),
    });
  }

  return out;
}

function sortedKey(a: string, b: string): string {
  return a < b ? `${a}/${b}` : `${b}/${a}`;
}

/** Upper-tail p-value of chi-square distribution: 1 - F(chi2, df). */
export function chiSquarePValue(chiSq: number, df: number): number {
  if (!Number.isFinite(chiSq) || chiSq < 0 || df <= 0) return NaN;
  return 1 - regularizedGammaP(df / 2, chiSq / 2);
}

/**
 * Regularized lower incomplete gamma function P(a, x) = γ(a,x) / Γ(a).
 * Series for x < a+1, continued fraction otherwise — Numerical Recipes 6.2.
 */
function regularizedGammaP(a: number, x: number): number {
  if (x < 0 || a <= 0) return NaN;
  if (x === 0) return 0;
  if (x < a + 1) return gammaSeries(a, x);
  return 1 - gammaContFrac(a, x);
}

function gammaSeries(a: number, x: number): number {
  const ITMAX = 200;
  const EPS = 3e-12;
  let ap = a;
  let sum = 1 / a;
  let del = sum;
  for (let n = 0; n < ITMAX; n++) {
    ap++;
    del *= x / ap;
    sum += del;
    if (Math.abs(del) < Math.abs(sum) * EPS) break;
  }
  return sum * Math.exp(-x + a * Math.log(x) - logGamma(a));
}

function gammaContFrac(a: number, x: number): number {
  const ITMAX = 200;
  const EPS = 3e-12;
  const FPMIN = 1e-300;
  let b = x + 1 - a;
  let c = 1 / FPMIN;
  let d = 1 / b;
  let h = d;
  for (let i = 1; i <= ITMAX; i++) {
    const an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = b + an / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1) < EPS) break;
  }
  return Math.exp(-x + a * Math.log(x) - logGamma(a)) * h;
}

/** Lanczos log-gamma. Accurate to ~14 digits. */
export function logGamma(z: number): number {
  const g = 7;
  const c = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7,
  ];
  if (z < 0.5) {
    return Math.log(Math.PI / Math.sin(Math.PI * z)) - logGamma(1 - z);
  }
  z -= 1;
  let x = c[0]!;
  for (let i = 1; i < g + 2; i++) x += c[i]! / (z + i);
  const t = z + g + 0.5;
  return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
}
