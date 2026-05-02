import type { Genotypes } from "../parser.ts";
import { distinctAlleles } from "../parser.ts";

/**
 * Allele frequency table for a population (or all individuals).
 *
 * Rows = alleles (union across loci, sorted numerically), columns = loci.
 * Entries are frequencies (allele count / total typed gene copies at locus),
 * or null if the allele was not observed at that locus. The trailing "N"
 * row gives the number of typed gene copies per locus.
 */
export interface FrequencyTable {
  loci: string[];
  /** Union of alleles across all loci, sorted. */
  alleles: string[];
  /** values[alleleIdx][locusIdx] = frequency or null */
  values: (number | null)[][];
  /** N[locusIdx] = number of typed gene copies (= 2*N_diploid - missing) */
  N: number[];
}

export function alleleFrequencies(
  genos: Genotypes,
  indMask?: boolean[],
): FrequencyTable {
  const counts: Map<string, number>[] = genos.loci.map(() => new Map());
  const N: number[] = genos.loci.map(() => 0);

  for (let l = 0; l < genos.loci.length; l++) {
    const tuples = genos.alleles[l]!;
    const c = counts[l]!;
    for (let i = 0; i < tuples.length; i++) {
      if (indMask && !indMask[i]) continue;
      for (const a of tuples[i]!) {
        if (a === null) continue;
        c.set(a, (c.get(a) ?? 0) + 1);
        N[l]!++;
      }
    }
  }

  // Union of alleles, sorted numerically when possible.
  const allAlleles = new Set<string>();
  for (const c of counts) for (const k of c.keys()) allAlleles.add(k);
  const alleles = Array.from(allAlleles).sort((a, b) => {
    const na = Number(a);
    const nb = Number(b);
    if (Number.isFinite(na) && Number.isFinite(nb)) return na - nb;
    return a.localeCompare(b);
  });

  const values: (number | null)[][] = alleles.map((a) =>
    counts.map((c, l) => {
      const cnt = c.get(a);
      if (cnt === undefined) return null;
      const n = N[l]!;
      return n === 0 ? null : cnt / n;
    }),
  );

  return { loci: genos.loci, alleles, values, N };
}

/** Per-locus allele frequencies as Maps (convenience for downstream stats). */
export function perLocusFrequencies(
  genos: Genotypes,
  indMask?: boolean[],
): { freqs: Map<string, number>[]; N: number[] } {
  const freqs: Map<string, number>[] = [];
  const N: number[] = [];
  for (let l = 0; l < genos.loci.length; l++) {
    const counts = new Map<string, number>();
    let n = 0;
    const tuples = genos.alleles[l]!;
    for (let i = 0; i < tuples.length; i++) {
      if (indMask && !indMask[i]) continue;
      for (const a of tuples[i]!) {
        if (a === null) continue;
        counts.set(a, (counts.get(a) ?? 0) + 1);
        n++;
      }
    }
    const f = new Map<string, number>();
    if (n > 0) for (const [k, v] of counts) f.set(k, v / n);
    freqs.push(f);
    N.push(n);
  }
  return { freqs, N };
}

/** Distinct alleles per locus (used by various stats). */
export function locusAlleles(genos: Genotypes): string[][] {
  return genos.loci.map((_, i) => distinctAlleles(genos, i));
}
