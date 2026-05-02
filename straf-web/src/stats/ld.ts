import type { Genotypes } from "../parser.ts";
import { chiSquarePValue } from "./hwe.ts";

/**
 * Pairwise linkage disequilibrium test between every pair of loci.
 *
 * Method: chi-square test of independence on the 2-locus *genotype* contingency
 * table, restricted to individuals fully typed at both loci.
 *
 *   H0: genotype at locus A is independent of genotype at locus B
 *   χ² = Σ_{ij} (O_ij − E_ij)² / E_ij,   E_ij = row_i · col_j / N
 *   df  = (R − 1)(C − 1)   where R, C are observed-genotype counts
 *
 * Notes:
 *   - This is an asymptotic approximation. The original STRAF uses Genepop's
 *     exact MCMC test (Slatkin 1994 / Raymond & Rousset 1995). The asymptotic
 *     test is conservative when many cells are sparse — be cautious with
 *     low-frequency genotypes.
 *   - p-values are clamped to (1/N, 1) where N is the per-pair sample size,
 *     mirroring what the original does after the MCMC test.
 *   - If the test is undefined (df ≤ 0 or N = 0), p is reported as NaN.
 */

export interface LdResult {
  locus1: string;
  locus2: string;
  N: number;
  df: number;
  chiSq: number;
  pValue: number;
}

export function pairwiseLd(
  genos: Genotypes,
  indMask?: boolean[],
): LdResult[] {
  const out: LdResult[] = [];
  for (let l1 = 0; l1 < genos.loci.length - 1; l1++) {
    for (let l2 = l1 + 1; l2 < genos.loci.length; l2++) {
      out.push(pairLd(genos, l1, l2, indMask));
    }
  }
  return out;
}

function genotypeKey(tuple: (string | null)[]): string | null {
  for (const a of tuple) if (a === null) return null;
  // Sort copies so unordered genotype.
  const sorted = [...(tuple as string[])].sort();
  return sorted.join("/");
}

function pairLd(
  genos: Genotypes,
  l1: number,
  l2: number,
  indMask?: boolean[],
): LdResult {
  const t1 = genos.alleles[l1]!;
  const t2 = genos.alleles[l2]!;
  // table[g1][g2] count
  const counts = new Map<string, Map<string, number>>();
  const rowTotal = new Map<string, number>();
  const colTotal = new Map<string, number>();
  let N = 0;
  for (let i = 0; i < t1.length; i++) {
    if (indMask && !indMask[i]) continue;
    const k1 = genotypeKey(t1[i]!);
    const k2 = genotypeKey(t2[i]!);
    if (k1 === null || k2 === null) continue;
    if (!counts.has(k1)) counts.set(k1, new Map());
    const row = counts.get(k1)!;
    row.set(k2, (row.get(k2) ?? 0) + 1);
    rowTotal.set(k1, (rowTotal.get(k1) ?? 0) + 1);
    colTotal.set(k2, (colTotal.get(k2) ?? 0) + 1);
    N++;
  }

  const R = rowTotal.size;
  const C = colTotal.size;
  const df = (R - 1) * (C - 1);

  if (N === 0 || df <= 0) {
    return {
      locus1: genos.loci[l1]!,
      locus2: genos.loci[l2]!,
      N,
      df,
      chiSq: NaN,
      pValue: NaN,
    };
  }

  let chi = 0;
  for (const [k1, row] of counts) {
    const r = rowTotal.get(k1)!;
    for (const [k2, o] of row) {
      const c = colTotal.get(k2)!;
      const e = (r * c) / N;
      if (e > 0) chi += ((o - e) ** 2) / e;
    }
    // Cells with O = 0 still contribute E²/E = E.
    for (const [k2, c] of colTotal) {
      if (!row.has(k2)) {
        const e = (r * c) / N;
        chi += e;
      }
    }
  }

  let p = chiSquarePValue(chi, df);
  if (Number.isFinite(p)) {
    p = Math.max(1 / Math.max(N, 1), Math.min(1, p));
  }

  return {
    locus1: genos.loci[l1]!,
    locus2: genos.loci[l2]!,
    N,
    df,
    chiSq: chi,
    pValue: p,
  };
}

/** Reshape pairwise LD results into a square symmetric p-value matrix. */
export function ldMatrix(
  results: LdResult[],
  loci: string[],
): { loci: string[]; matrix: (number | null)[][] } {
  const idx = new Map(loci.map((l, i) => [l, i]));
  const matrix: (number | null)[][] = loci.map(() => loci.map(() => null));
  for (let i = 0; i < loci.length; i++) matrix[i]![i] = NaN;
  for (const r of results) {
    const i = idx.get(r.locus1)!;
    const j = idx.get(r.locus2)!;
    matrix[i]![j] = r.pValue;
    matrix[j]![i] = r.pValue;
  }
  return { loci, matrix };
}
