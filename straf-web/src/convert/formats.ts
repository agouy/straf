/**
 * STRAF file conversion utilities.
 *
 * Each function takes the parsed genotype data and produces a string in the
 * target format. Output strings match the originals byte-for-byte where
 * practical, accepting differences only in trailing whitespace and date/title
 * stamps.
 *
 * Supported targets:
 *   - Genepop  (microsatellite, 3-digit allele encoding)
 *   - Arlequin (microsatellite, MICROSAT data type)
 *   - Familias (per-locus allele frequency table)
 *   - Euroformix (CSV allele frequency table without the trailing N row)
 *   - STRmix   (CSV allele frequency table including the N row)
 *   - LRmix    (CSV allele frequency table, NA encoded as 0)
 */

import type { Genotypes } from "../parser.ts";
import { alleleFrequencies } from "../stats/freq.ts";
import { popMask, allMask } from "../parser.ts";

/**
 * Convert a numeric STR allele label to the 3-digit zero-padded encoding used
 * by Genepop, Arlequin, and friends. Allele "9.3" → "093" (drop the dot,
 * pad to 3 digits). Plain integers like "12" → "120" (3 digits with the
 * fractional decade as 0).
 *
 * Mirrors the ad-hoc transform inside straf2genepop / straf2arlequin in the
 * original R code.
 */
export function encodeAllele3Digit(a: string): string {
  const compact = a.replace(".", "");
  if (compact.length === 1) return compact + "00";
  if (compact.length === 2) return compact + "0";
  if (compact.length === 3) return compact;
  throw new Error(`Cannot encode allele "${a}" in the 3-digit format.`);
}

// --- Genepop ---------------------------------------------------------------

export function toGenepop(genos: Genotypes): string {
  const lines: string[] = ["STRAF-generated GENEPOP input file."];
  for (const l of genos.loci) lines.push(l);

  // Group rows by population, preserving first-seen order.
  const popOrder: string[] = [];
  const byPop: Map<string, number[]> = new Map();
  for (let i = 0; i < genos.populations.length; i++) {
    const p = genos.populations[i]!;
    if (!byPop.has(p)) {
      popOrder.push(p);
      byPop.set(p, []);
    }
    byPop.get(p)!.push(i);
  }

  for (const p of popOrder) {
    lines.push("Pop");
    for (const idx of byPop.get(p)!) {
      const ind = genos.individuals[idx]!;
      const alleles: string[] = [];
      for (let l = 0; l < genos.loci.length; l++) {
        const tuple = genos.alleles[l]![idx]!;
        // Genepop expects ploidy concatenated alleles per locus.
        const encoded = tuple.map((a) => (a === null ? "000" : encodeAllele3Digit(a)));
        alleles.push(encoded.join(""));
      }
      lines.push(`${ind}\t,\t${alleles.join("\t")}`);
    }
  }

  return lines.join("\n") + "\n";
}

// --- Arlequin --------------------------------------------------------------

export function toArlequin(genos: Genotypes): string {
  if (genos.ploidy !== 2) {
    throw new Error("Arlequin export currently supports diploid data only.");
  }

  const popOrder: string[] = [];
  const byPop: Map<string, number[]> = new Map();
  for (let i = 0; i < genos.populations.length; i++) {
    const p = genos.populations[i]!;
    if (!byPop.has(p)) {
      popOrder.push(p);
      byPop.set(p, []);
    }
    byPop.get(p)!.push(i);
  }

  const out: string[] = [];
  out.push("[Profile]");
  out.push('Title="STRAF-generated Arlequin file."');
  out.push(`NbSamples=${popOrder.length}`);
  out.push("DataType=MICROSAT");
  out.push("");
  out.push("GenotypicData=1");
  out.push("GameticPhase=0");
  out.push('MissingData="?"');
  out.push("LocusSeparator=WHITESPACE");
  out.push("");
  out.push("[Data]");
  out.push("");
  out.push("[[Samples]]");
  out.push("");

  for (const p of popOrder) {
    const idxs = byPop.get(p)!;
    out.push(`SampleName="${p}"`);
    out.push(`SampleSize=${idxs.length}`);
    out.push("SampleData={");
    for (const i of idxs) {
      const ind = genos.individuals[i]!;
      const a1: string[] = [];
      const a2: string[] = [];
      for (let l = 0; l < genos.loci.length; l++) {
        const tuple = genos.alleles[l]![i]!;
        a1.push(tuple[0] === null ? "?" : encodeAllele3Digit(tuple[0]!));
        a2.push(tuple[1] === null ? "?" : encodeAllele3Digit(tuple[1]!));
      }
      out.push([ind, "1", ...a1].join("\t"));
      out.push(["", "", ...a2].join("\t"));
    }
    out.push("}");
    out.push("");
  }

  return out.join("\n") + "\n";
}

// --- Familias --------------------------------------------------------------

/**
 * Familias allele-frequency-database format: one block per locus, lines of
 * "<allele>\t<frequency>". Blocks separated by a blank line.
 */
export function toFamilias(genos: Genotypes, pop = "all"): string {
  const mask = pop === "all" ? allMask(genos) : popMask(genos, pop);
  const tbl = alleleFrequencies(genos, mask);
  const blocks: string[] = [];
  for (let l = 0; l < tbl.loci.length; l++) {
    const lines: string[] = [tbl.loci[l]!];
    for (let a = 0; a < tbl.alleles.length; a++) {
      const f = tbl.values[a]![l];
      if (f === null || f === undefined) continue;
      lines.push(`${tbl.alleles[a]}\t${f}`);
    }
    blocks.push(lines.join("\n"));
  }
  return blocks.join("\n\n") + "\n";
}

// --- Allele-frequency CSV exports (Euroformix / STRmix / LRmix) ------------

/**
 * Allele × locus frequency table as CSV. Variants:
 *   - Euroformix: drop the trailing N row, encode missing as empty string.
 *   - STRmix:    keep the N row, encode missing as empty string.
 *   - LRmix:     drop the N row, encode missing as "0".
 */
export type CsvFreqVariant = "euroformix" | "strmix" | "lrmix";

export function toAlleleFreqCsv(
  genos: Genotypes,
  variant: CsvFreqVariant,
  pop = "all",
): string {
  const mask = pop === "all" ? allMask(genos) : popMask(genos, pop);
  const tbl = alleleFrequencies(genos, mask);

  const missing = variant === "lrmix" ? "0" : "";
  const dropN = variant === "euroformix" || variant === "lrmix";

  const headers = ["Allele", ...tbl.loci];
  const rows: string[][] = tbl.alleles.map((a, ai) => [
    a,
    ...tbl.values[ai]!.map((v) => (v === null || v === undefined ? missing : String(v))),
  ]);
  if (!dropN) rows.push(["N", ...tbl.N.map((n) => String(n))]);

  return [headers.join(","), ...rows.map((r) => r.join(","))].join("\n") + "\n";
}
