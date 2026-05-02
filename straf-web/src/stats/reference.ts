/**
 * Parse the STRidER-style reference frequency CSV used by the original STRAF.
 *
 * Format: blocks of one locus each, separated by a single-field line giving
 * the locus name. Within a block:
 *
 *   D1S1656                                           (locus name, 1 field)
 *   Allele,"Entire Database",Africa,...               (header: pop names)
 *   ,13966,1243,...                                   (sample sizes; allele blank)
 *   1,0.000967,,,...                                  (allele rows)
 *   8,0.000358,0.003218,,...
 *   ...
 *
 * Cells separated by commas. Quoted fields (with embedded commas) are
 * supported. Empty cells = missing.
 *
 * The output is a population × (locus_allele) frequency matrix, with column
 * names of the form `<locus>_<allele>` and rows labeled by population.
 * Missing cells are stored as null.
 */

export interface ReferenceFrequencies {
  populations: string[];
  /** columns[k] = "LOCUS_ALLELE" string */
  columns: string[];
  /** values[popIdx][colIdx] = frequency or null */
  values: (number | null)[][];
}

export function parseReferenceCsv(text: string): ReferenceFrequencies {
  const lines = text
    .split(/\r?\n/)
    .map((l) => l.replace(/^﻿/, ""))
    .filter((l) => l.length > 0);

  // Identify locus blocks: a "name line" has exactly one non-empty field.
  const blocks: { name: string; start: number; end: number }[] = [];
  for (let i = 0; i < lines.length; i++) {
    const cells = parseCsvLine(lines[i]!);
    const nonEmpty = cells.filter((c) => c.trim().length > 0);
    if (nonEmpty.length === 1 && cells.length === 1) {
      // Possibly a locus name. Confirm next line starts with "Allele".
      if (i + 1 < lines.length) {
        const next = parseCsvLine(lines[i + 1]!);
        if (next[0]?.trim().toLowerCase() === "allele") {
          if (blocks.length > 0) blocks[blocks.length - 1]!.end = i;
          blocks.push({ name: nonEmpty[0]!, start: i, end: lines.length });
        }
      }
    }
  }

  if (blocks.length === 0) {
    throw new Error("Reference CSV: no locus blocks recognised.");
  }

  // Collect populations from the first block's header.
  const firstHeader = parseCsvLine(lines[blocks[0]!.start + 1]!);
  const populations = firstHeader.slice(1).map((s) => s.trim());

  const columns: string[] = [];
  const valuesByCol: (number | null)[][] = []; // valuesByCol[colIdx][popIdx]

  for (const block of blocks) {
    const header = parseCsvLine(lines[block.start + 1]!);
    const pops = header.slice(1).map((s) => s.trim());
    // Build mapping from this block's populations to the canonical order.
    // (In practice all blocks in the STRidER file share the same population order, but we don't assume it.)
    const popMap: number[] = pops.map((p) => populations.indexOf(p));
    // Sample-size line is start + 2 (allele cell empty). Skip it.
    let dataStart = block.start + 2;
    const sizeRow = parseCsvLine(lines[dataStart]!);
    if (sizeRow[0]?.trim() !== "") {
      // No sample-size row — first data row starts at +2 anyway.
      dataStart = block.start + 2;
    } else {
      dataStart = block.start + 3;
    }

    for (let r = dataStart; r < block.end; r++) {
      const cells = parseCsvLine(lines[r]!);
      const allele = cells[0]?.trim();
      if (!allele) continue;
      const colName = `${block.name}_${allele}`;
      const col: (number | null)[] = populations.map(() => null);
      for (let k = 0; k < pops.length; k++) {
        const target = popMap[k];
        if (target === undefined || target < 0) continue;
        const raw = cells[k + 1]?.trim() ?? "";
        if (raw === "") continue;
        const v = Number(raw);
        if (Number.isFinite(v)) col[target] = v;
      }
      columns.push(colName);
      valuesByCol.push(col);
    }
  }

  // Transpose to rows = populations.
  const values: (number | null)[][] = populations.map((_, p) =>
    valuesByCol.map((col) => col[p]!),
  );

  return { populations, columns, values };
}

/**
 * Drop columns where any selected population has null. Mirrors the
 * `colSums(is.na(X)) == 0` filter in the original R code.
 */
export function dropMissingColumns(
  ref: ReferenceFrequencies,
  selectedPops?: string[],
): ReferenceFrequencies {
  const keepRow = ref.populations.map((p) => !selectedPops || selectedPops.includes(p));
  const keep: boolean[] = ref.columns.map((_, c) => {
    for (let p = 0; p < ref.populations.length; p++) {
      if (!keepRow[p]) continue;
      if (ref.values[p]![c] === null) return false;
    }
    return true;
  });
  const newCols = ref.columns.filter((_, c) => keep[c]);
  const newValues = ref.values.map((row) => row.filter((_, c) => keep[c]));
  return { populations: ref.populations, columns: newCols, values: newValues };
}

/**
 * Nei-style distance from a population × allele-frequency matrix.
 *
 *   I_ij = Σ_a p_ia · p_ja  /  sqrt(Σ_a p_ia² · Σ_a p_ja²)
 *   D_ij = -log(I_ij)
 *
 * (This matches the formula `d <- X %*% t(X) / sqrt(diag · diag); -log(d)`
 *  used by STRAF on the reference table.)
 */
export function neiDistanceFromMatrix(values: number[][]): number[][] {
  const n = values.length;
  const matrix: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
  const dot = (a: number[], b: number[]) => {
    let s = 0;
    for (let k = 0; k < a.length; k++) s += a[k]! * b[k]!;
    return s;
  };
  const norms = values.map((row) => Math.sqrt(dot(row, row)));
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const num = dot(values[i]!, values[j]!);
      const denom = norms[i]! * norms[j]!;
      const I = denom === 0 ? 0 : num / denom;
      const d = -Math.log(Math.max(1e-12, Math.min(1, I)));
      matrix[i]![j] = d;
      matrix[j]![i] = d;
    }
  }
  return matrix;
}

// --- CSV line parser (RFC-4180-ish, handles quoted fields and embedded quotes) ---
function parseCsvLine(line: string): string[] {
  const out: string[] = [];
  let cur = "";
  let inQuotes = false;
  for (let i = 0; i < line.length; i++) {
    const ch = line[i]!;
    if (inQuotes) {
      if (ch === '"') {
        if (line[i + 1] === '"') {
          cur += '"';
          i++;
        } else {
          inQuotes = false;
        }
      } else {
        cur += ch;
      }
    } else {
      if (ch === ",") {
        out.push(cur);
        cur = "";
      } else if (ch === '"' && cur === "") {
        inQuotes = true;
      } else {
        cur += ch;
      }
    }
  }
  out.push(cur);
  return out;
}
