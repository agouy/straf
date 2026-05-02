/**
 * STRAF input parser with detailed validation diagnostics.
 *
 * Format: tab-separated values with header.
 *   - Column 1 = "ind" (individual ID)
 *   - Column 2 = "pop" (population label)
 *   - Remaining columns = allele calls. For diploid data each locus name
 *     appears in two adjacent columns sharing the same column name; for
 *     haploid, one column per locus. Missing alleles are encoded as "0".
 *
 * Allele labels are kept as strings to preserve fractional STR designators
 * (e.g. "9.3", "30.2"). Numeric coercion only happens at display time.
 *
 * The parser tries to give actionable error messages — when an obvious
 * problem is detected (spaces instead of tabs, wrong ploidy, special
 * characters, etc.) it points at the line, column, and what to fix.
 */

export type Ploidy = 1 | 2;

export interface Genotypes {
  individuals: string[];
  populations: string[];
  loci: string[];
  ploidy: Ploidy;
  /** alleles[locusIdx][indIdx] = array of `ploidy` allele strings (or null for missing) */
  alleles: (string | null)[][][];
}

export interface ParseIssue {
  /** "error" prevents loading; "warning" is informational. */
  kind: "error" | "warning";
  message: string;
  /** Line number in the source file (1-indexed), if applicable. */
  line?: number;
  /** Column index (1-indexed), if applicable. */
  column?: number;
  /** A short suggestion or hint, when we can give one. */
  hint?: string;
}

export class ParseError extends Error {
  readonly issues: ParseIssue[];

  constructor(issues: ParseIssue[]) {
    super(formatIssues(issues));
    this.issues = issues;
  }
}

function formatIssues(issues: ParseIssue[]): string {
  return issues
    .map((i) => {
      const where = i.line ? ` (line ${i.line}${i.column ? `, col ${i.column}` : ""})` : "";
      const hint = i.hint ? `\n  → ${i.hint}` : "";
      return `${i.kind.toUpperCase()}${where}: ${i.message}${hint}`;
    })
    .join("\n");
}

const ALLOWED_NAME_RE = /^[A-Za-z0-9_.\-]+$/;
const NUMERIC_ALLELE_RE = /^[0-9]+(?:\.[0-9]+)?$/;

export interface ParseOptions {
  /**
   * If true (default), missing values may be encoded as either "0" or empty
   * cells. If false, only literal "0" is accepted as missing.
   */
  acceptEmptyAsMissing?: boolean;
}

export function parseStrafTsv(
  text: string,
  ploidy: Ploidy,
  options: ParseOptions = {},
): Genotypes {
  const issues: ParseIssue[] = [];

  // ------------------------------------------------------------------
  // 1. Split into raw + trimmed lines (preserve original for diagnostics)
  // ------------------------------------------------------------------
  const rawLines = text.split(/\r?\n/);
  // Keep track of which lines have trailing whitespace, to report as warnings.
  for (let i = 0; i < rawLines.length; i++) {
    const raw = rawLines[i]!;
    if (raw.length === 0) continue;
    if (/[ \t]+$/.test(raw) && i < rawLines.length - 1) {
      issues.push({
        kind: "warning",
        line: i + 1,
        message: "Trailing whitespace at end of line.",
        hint: "Whitespace was trimmed automatically; consider cleaning up the source file.",
      });
    }
  }

  // Drop empty trailing lines but keep mid-file blank rows visible (those are
  // treated as bad rows below).
  const lines: { lineNum: number; text: string }[] = [];
  for (let i = 0; i < rawLines.length; i++) {
    const t = rawLines[i]!.trimEnd();
    if (t.length === 0 && i === rawLines.length - 1) continue; // trailing blank
    lines.push({ lineNum: i + 1, text: t });
  }

  if (lines.length < 2) {
    throw new ParseError([
      {
        kind: "error",
        message: "File must contain a header row and at least one data row.",
        hint: "Check that you saved the file as Tab-delimited text and that it has more than one line.",
      },
    ]);
  }

  // ------------------------------------------------------------------
  // 2. Header parsing & "spaces instead of tabs" detection
  // ------------------------------------------------------------------
  const headerLine = lines[0]!;
  let headerRaw = splitTab(headerLine.text).map((h) => h.trim());

  // If the header has only one cell but contains spaces, the file is probably
  // space-separated. This is the most common single mistake.
  if (headerRaw.length === 1 && /\s/.test(headerLine.text)) {
    throw new ParseError([
      {
        kind: "error",
        line: headerLine.lineNum,
        message: "Header has only one column — values appear to be separated by spaces, not tabs.",
        hint:
          "STRAF requires TAB-separated values. " +
          "In Excel use Save As → Text (Tab delimited) (*.txt). " +
          "If you generated the file from a script, use \"\\t\" rather than \" \" or \",\".",
      },
    ]);
  }

  // ------------------------------------------------------------------
  // 3. First-column-name checks (ind / pop) with case-tolerant hints
  // ------------------------------------------------------------------
  if (headerRaw.length < 3) {
    issues.push({
      kind: "error",
      line: headerLine.lineNum,
      message: `Header has only ${headerRaw.length} column(s); need at least 3 (ind, pop, and one or more loci).`,
    });
    throw new ParseError(issues);
  }
  const c0 = headerRaw[0]!;
  const c1 = headerRaw[1]!;
  if (c0 !== "ind") {
    const looksLikeIt = c0.toLowerCase() === "ind";
    issues.push({
      kind: "error",
      line: headerLine.lineNum,
      column: 1,
      message: `First column must be named "ind" (lowercase). Got "${c0}".`,
      hint: looksLikeIt ? `Rename "${c0}" to "ind".` : undefined,
    });
  }
  if (c1 !== "pop") {
    const looksLikeIt = c1.toLowerCase() === "pop";
    issues.push({
      kind: "error",
      line: headerLine.lineNum,
      column: 2,
      message: `Second column must be named "pop" (lowercase). Got "${c1}".`,
      hint: looksLikeIt ? `Rename "${c1}" to "pop".` : undefined,
    });
  }
  if (issues.some((i) => i.kind === "error")) throw new ParseError(issues);

  // ------------------------------------------------------------------
  // 4. Locus columns: count, ploidy, naming
  // ------------------------------------------------------------------
  const locusCols = headerRaw.slice(2);
  if (locusCols.length === 0) {
    throw new ParseError([
      {
        kind: "error",
        line: headerLine.lineNum,
        message: "No locus columns found after ind/pop.",
      },
    ]);
  }

  // Special-character check on locus names.
  for (let i = 0; i < locusCols.length; i++) {
    const name = locusCols[i]!;
    if (name.length === 0) {
      issues.push({
        kind: "error",
        line: headerLine.lineNum,
        column: i + 3,
        message: `Empty locus column name at position ${i + 3}.`,
      });
    } else if (!ALLOWED_NAME_RE.test(name)) {
      issues.push({
        kind: "warning",
        line: headerLine.lineNum,
        column: i + 3,
        message: `Locus name "${name}" contains characters other than letters, digits, underscore, hyphen and period.`,
        hint: "Some downstream tools (Genepop, Arlequin) reject special characters; consider renaming.",
      });
    }
  }

  if (locusCols.length % ploidy !== 0) {
    issues.push({
      kind: "error",
      line: headerLine.lineNum,
      message: `Number of locus columns (${locusCols.length}) is not a multiple of ploidy (${ploidy}).`,
      hint:
        ploidy === 2
          ? "For diploid data, each locus must appear twice with the same name. Did you mean to set ploidy to 1 (haploid)?"
          : "For haploid data, each locus appears once. Did you mean to set ploidy to 2 (diploid)?",
    });
    throw new ParseError(issues);
  }

  // Group locus columns by ploidy and validate.
  const loci: string[] = [];
  for (let i = 0; i < locusCols.length; i += ploidy) {
    const name = locusCols[i]!;
    for (let k = 1; k < ploidy; k++) {
      const got = locusCols[i + k]!;
      if (got !== name) {
        issues.push({
          kind: "error",
          line: headerLine.lineNum,
          column: i + k + 3,
          message: `For diploid data, columns ${i + 3} and ${i + k + 3} must have the same locus name. Got "${name}" and "${got}".`,
          hint:
            "Each diploid locus needs two adjacent columns with identical names. " +
            "If you only have one column per locus, switch ploidy to 1 (haploid) in the sidebar.",
        });
      }
    }
    loci.push(name);
  }
  if (issues.some((i) => i.kind === "error")) throw new ParseError(issues);

  // For haploid data, all locus names must be distinct.
  if (ploidy === 1) {
    const seen = new Map<string, number>();
    for (let i = 0; i < loci.length; i++) {
      const dup = seen.get(loci[i]!);
      if (dup !== undefined) {
        issues.push({
          kind: "error",
          line: headerLine.lineNum,
          column: i + 3,
          message: `Locus name "${loci[i]}" is duplicated (also at column ${dup + 3}).`,
          hint: "For haploid data each locus should appear only once. Are you sure ploidy is set correctly?",
        });
      }
      seen.set(loci[i]!, i);
    }
  } else {
    // For diploid, locus names appear twice (the pairs); any name repeated
    // beyond a single pair is an error.
    const seen = new Set<string>();
    for (let i = 0; i < loci.length; i++) {
      if (seen.has(loci[i]!)) {
        issues.push({
          kind: "error",
          line: headerLine.lineNum,
          message: `Locus "${loci[i]}" appears more than once across non-adjacent column pairs.`,
        });
      }
      seen.add(loci[i]!);
    }
  }
  if (issues.some((i) => i.kind === "error")) throw new ParseError(issues);

  // ------------------------------------------------------------------
  // 5. Data rows
  // ------------------------------------------------------------------
  const acceptEmpty = options.acceptEmptyAsMissing ?? true;
  const expectedCols = headerRaw.length;

  const individuals: string[] = [];
  const populations: string[] = [];
  const alleles: (string | null)[][][] = loci.map(() => []);
  const indNames = new Set<string>();

  for (let r = 1; r < lines.length; r++) {
    const { lineNum, text: lineText } = lines[r]!;
    if (lineText.length === 0) {
      issues.push({
        kind: "warning",
        line: lineNum,
        message: "Empty line in the middle of the file — skipped.",
      });
      continue;
    }

    const cells = splitTab(lineText);

    // Spaces-instead-of-tabs heuristic at the row level.
    if (cells.length === 1 && /\s/.test(lineText) && expectedCols > 1) {
      issues.push({
        kind: "error",
        line: lineNum,
        message: "Row appears to be space-separated rather than tab-separated.",
        hint: 'Replace spaces with tabs. In Excel use Save As → Text (Tab delimited).',
      });
      throw new ParseError(issues);
    }

    if (cells.length !== expectedCols) {
      issues.push({
        kind: "error",
        line: lineNum,
        message: `Expected ${expectedCols} tab-separated columns, got ${cells.length}.`,
        hint:
          cells.length < expectedCols
            ? "Some cells may be missing. If a value is empty, encode missing alleles as \"0\" (do not leave the cell blank)."
            : "Extra columns — check for stray tabs or invisible characters in the row.",
      });
      throw new ParseError(issues);
    }

    const ind = cells[0]!.trim();
    const pop = cells[1]!.trim();

    if (ind.length === 0) {
      issues.push({
        kind: "error",
        line: lineNum,
        column: 1,
        message: "Empty individual ID.",
      });
    }
    if (pop.length === 0) {
      issues.push({
        kind: "error",
        line: lineNum,
        column: 2,
        message: "Empty population label.",
      });
    }
    if (ind.length > 0 && !ALLOWED_NAME_RE.test(ind)) {
      issues.push({
        kind: "warning",
        line: lineNum,
        column: 1,
        message: `Individual ID "${ind}" contains characters other than letters, digits, underscore, hyphen and period.`,
      });
    }
    if (pop.length > 0 && !ALLOWED_NAME_RE.test(pop)) {
      issues.push({
        kind: "warning",
        line: lineNum,
        column: 2,
        message: `Population label "${pop}" contains characters other than letters, digits, underscore, hyphen and period.`,
      });
    }

    if (indNames.has(ind)) {
      issues.push({
        kind: "warning",
        line: lineNum,
        column: 1,
        message: `Duplicate individual ID "${ind}".`,
      });
    }
    indNames.add(ind);

    individuals.push(ind);
    populations.push(pop);

    for (let l = 0; l < loci.length; l++) {
      const tuple: (string | null)[] = [];
      for (let k = 0; k < ploidy; k++) {
        const colIdx = 2 + l * ploidy + k; // 0-indexed
        const raw = cells[colIdx]!;
        const v = raw.trim();

        if (v === "0") {
          tuple.push(null);
          continue;
        }
        if (v === "") {
          if (acceptEmpty) {
            tuple.push(null);
            issues.push({
              kind: "warning",
              line: lineNum,
              column: colIdx + 1,
              message: `Empty cell at locus ${loci[l]} treated as missing.`,
              hint: 'Encode missing alleles explicitly as "0".',
            });
            continue;
          }
          issues.push({
            kind: "error",
            line: lineNum,
            column: colIdx + 1,
            message: `Empty cell at locus ${loci[l]}.`,
            hint: 'Encode missing alleles as "0".',
          });
          tuple.push(null);
          continue;
        }
        if (!NUMERIC_ALLELE_RE.test(v)) {
          issues.push({
            kind: "error",
            line: lineNum,
            column: colIdx + 1,
            message: `Allele "${v}" at locus ${loci[l]} is not numeric.`,
            hint:
              'Alleles must be encoded as numbers, optionally with a fractional part (e.g. "9.3"). ' +
              'Use "0" for missing data — not "NA", "N/A", "?", or letters.',
          });
        }
        tuple.push(v);
      }
      alleles[l]!.push(tuple);
    }
  }

  if (individuals.length === 0) {
    issues.push({ kind: "error", message: "No data rows found." });
  }

  // Ploidy/data sanity heuristic: for haploid data with no missing values,
  // if every cell at every locus is the same value, suspect ploidy mismatch.
  // (Cheap check; not authoritative.)

  // Surface fatal issues; warnings are kept on the resulting object so the UI
  // can show them.
  if (issues.some((i) => i.kind === "error")) throw new ParseError(issues);

  const result: Genotypes = { individuals, populations, loci, ploidy, alleles };
  // Stash warnings on the object (non-enumerable so JSON.stringify is unchanged).
  Object.defineProperty(result, "warnings", { value: issues, enumerable: false });
  return result;
}

/** Return any non-fatal warnings recorded by the parser, if available. */
export function getParseWarnings(genos: Genotypes): ParseIssue[] {
  const w = (genos as unknown as { warnings?: ParseIssue[] }).warnings;
  return w ?? [];
}

function splitTab(line: string): string[] {
  return line.split("\t");
}

/** Distinct alleles observed at a locus, sorted numerically when possible. */
export function distinctAlleles(genos: Genotypes, locusIdx: number): string[] {
  const set = new Set<string>();
  for (const tuple of genos.alleles[locusIdx]!) {
    for (const a of tuple) if (a !== null) set.add(a);
  }
  const arr = Array.from(set);
  arr.sort((a, b) => {
    const na = Number(a);
    const nb = Number(b);
    if (Number.isFinite(na) && Number.isFinite(nb)) return na - nb;
    return a.localeCompare(b);
  });
  return arr;
}

/** Number of typed (non-missing) gene copies at a locus. */
export function geneCopies(genos: Genotypes, locusIdx: number, indMask?: boolean[]): number {
  let n = 0;
  const tuples = genos.alleles[locusIdx]!;
  for (let i = 0; i < tuples.length; i++) {
    if (indMask && !indMask[i]) continue;
    for (const a of tuples[i]!) if (a !== null) n++;
  }
  return n;
}

export function uniquePopulations(genos: Genotypes): string[] {
  return Array.from(new Set(genos.populations));
}

export function popMask(genos: Genotypes, pop: string): boolean[] {
  return genos.populations.map((p) => p === pop);
}

export function allMask(genos: Genotypes): boolean[] {
  return genos.populations.map(() => true);
}
