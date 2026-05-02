import { test } from "node:test";
import { strict as assert } from "node:assert";

import {
  parseStrafTsv,
  ParseError,
  getParseWarnings,
} from "./parser.ts";
import type { ParseIssue } from "./parser.ts";

// --- helpers ---------------------------------------------------------------

function expectParseError(text: string, ploidy: 1 | 2): ParseIssue[] {
  try {
    parseStrafTsv(text, ploidy);
  } catch (e) {
    if (e instanceof ParseError) return e.issues;
    throw e;
  }
  throw new Error("expected ParseError, got success");
}

function lines(...rows: string[]): string {
  return rows.join("\n");
}

// --- happy paths -----------------------------------------------------------

test("parses minimal diploid data", () => {
  const txt = lines("ind\tpop\tL1\tL1", "i1\tA\t10\t11", "i2\tA\t10\t10");
  const g = parseStrafTsv(txt, 2);
  assert.equal(g.individuals.length, 2);
  assert.equal(g.loci.length, 1);
  assert.deepEqual(g.alleles[0], [["10", "11"], ["10", "10"]]);
});

test("parses minimal haploid data", () => {
  const txt = lines("ind\tpop\tL1\tL2", "i1\tA\t10\t12", "i2\tA\t11\t13");
  const g = parseStrafTsv(txt, 1);
  assert.equal(g.loci.length, 2);
  assert.deepEqual(g.alleles[0], [["10"], ["11"]]);
});

test("treats both 0 and empty as missing (with warning for empty)", () => {
  const txt = lines("ind\tpop\tL1\tL1", "i1\tA\t10\t0", "i2\tA\t\t11");
  const g = parseStrafTsv(txt, 2);
  assert.deepEqual(g.alleles[0]![0], ["10", null]);
  assert.deepEqual(g.alleles[0]![1], [null, "11"]);
  const w = getParseWarnings(g);
  assert.ok(w.some((i) => /empty cell/i.test(i.message)));
});

test("preserves point alleles like 9.3", () => {
  const txt = lines("ind\tpop\tFGA\tFGA", "i1\tA\t9.3\t11");
  const g = parseStrafTsv(txt, 2);
  assert.equal(g.alleles[0]![0]![0], "9.3");
});

test("trailing newline is OK", () => {
  const txt = "ind\tpop\tL1\tL1\ni1\tA\t10\t11\n";
  const g = parseStrafTsv(txt, 2);
  assert.equal(g.individuals.length, 1);
});

// --- 1. ind/pop column names ----------------------------------------------

test("error: missing ind / pop header names", () => {
  const txt = lines("Ind\tPOP\tL1\tL1", "i1\tA\t10\t11");
  const issues = expectParseError(txt, 2);
  assert.ok(issues.some((i) => /first column.*"ind"/i.test(i.message)));
  assert.ok(issues.some((i) => /second column.*"pop"/i.test(i.message)));
  // Hint should suggest renaming since the case difference looks intentional.
  assert.ok(issues.some((i) => /Rename "Ind" to "ind"/.test(i.hint ?? "")));
});

test("error: only 2 columns (no loci)", () => {
  const txt = lines("ind\tpop", "i1\tA");
  const issues = expectParseError(txt, 2);
  assert.ok(issues.some((i) => /at least 3/i.test(i.message)));
});

// --- 2. Locus naming for ploidy ------------------------------------------

test("error: diploid locus columns must repeat", () => {
  const txt = lines("ind\tpop\tL1\tL2", "i1\tA\t10\t11");
  const issues = expectParseError(txt, 2);
  assert.ok(
    issues.some((i) => /must have the same locus name/.test(i.message)),
    `got: ${JSON.stringify(issues)}`,
  );
  // Should hint at flipping ploidy
  assert.ok(issues.some((i) => /haploid/i.test(i.hint ?? "")));
});

test("error: haploid locus names must be distinct", () => {
  const txt = lines("ind\tpop\tL1\tL1", "i1\tA\t10\t11");
  const issues = expectParseError(txt, 1);
  assert.ok(issues.some((i) => /duplicated/i.test(i.message)));
});

test("error: odd column count for diploid", () => {
  const txt = lines("ind\tpop\tL1\tL1\tL2", "i1\tA\t10\t11\t12");
  const issues = expectParseError(txt, 2);
  assert.ok(issues.some((i) => /multiple of ploidy/i.test(i.message)));
});

// --- 3. Missing data encoding ---------------------------------------------

test('non-numeric allele "NA" is rejected with helpful hint', () => {
  const txt = lines("ind\tpop\tL1\tL1", "i1\tA\t10\tNA");
  const issues = expectParseError(txt, 2);
  assert.ok(issues.some((i) => /not numeric/i.test(i.message)));
  assert.ok(issues.some((i) => /\"NA\"|N\/A|0.+missing/i.test(i.hint ?? "")));
});

test('non-numeric allele "?" is rejected', () => {
  const txt = lines("ind\tpop\tL1\tL1", "i1\tA\t10\t?");
  const issues = expectParseError(txt, 2);
  assert.ok(issues.some((i) => /not numeric/i.test(i.message)));
});

// --- 4. Special characters in names ---------------------------------------

test("warning: special chars in locus name", () => {
  const txt = lines("ind\tpop\tD3*S1358\tD3*S1358", "i1\tA\t10\t11");
  const g = parseStrafTsv(txt, 2);
  const w = getParseWarnings(g);
  assert.ok(w.some((i) => /special|other than/i.test(i.message)));
});

test("warning: special chars in individual ID", () => {
  const txt = lines("ind\tpop\tL1\tL1", "i#1\tA\t10\t11");
  const g = parseStrafTsv(txt, 2);
  const w = getParseWarnings(g);
  assert.ok(w.some((i) => /Individual ID/i.test(i.message)));
});

// --- 5. Trailing whitespace ----------------------------------------------

test("warning: trailing whitespace on data row", () => {
  // Note: trimEnd happens before split; we want to see the warning *and*
  // still parse cleanly.
  const txt = "ind\tpop\tL1\tL1\ni1\tA\t10\t11   \ni2\tA\t11\t11";
  const g = parseStrafTsv(txt, 2);
  const w = getParseWarnings(g);
  assert.ok(w.some((i) => /trailing whitespace/i.test(i.message)));
  assert.equal(g.individuals.length, 2);
});

// --- 6. Spaces vs tabs ----------------------------------------------------

test("error: header is space-separated, not tabs", () => {
  const txt = "ind pop L1 L1\ni1 A 10 11";
  const issues = expectParseError(txt, 2);
  assert.ok(
    issues.some((i) => /space|tab/i.test(i.message)),
    `got: ${JSON.stringify(issues)}`,
  );
  assert.ok(issues.some((i) => /Excel|Tab delimited/i.test(i.hint ?? "")));
});

test("error: row has wrong column count", () => {
  const txt = lines("ind\tpop\tL1\tL1", "i1\tA\t10");
  const issues = expectParseError(txt, 2);
  assert.ok(issues.some((i) => /Expected 4 tab-separated/i.test(i.message)));
  assert.ok(issues[0]!.line === 2);
});

// --- 7. Empty / malformed lines ------------------------------------------

test("warning: empty mid-file row is skipped", () => {
  const txt = lines("ind\tpop\tL1\tL1", "i1\tA\t10\t11", "", "i2\tA\t11\t12");
  const g = parseStrafTsv(txt, 2);
  assert.equal(g.individuals.length, 2);
  const w = getParseWarnings(g);
  assert.ok(w.some((i) => /empty line/i.test(i.message)));
});

test("error: file with only a header", () => {
  const txt = "ind\tpop\tL1\tL1";
  const issues = expectParseError(txt, 2);
  assert.ok(issues.some((i) => /header.*and at least one data row/i.test(i.message)));
});

test("warning: duplicate individual IDs", () => {
  const txt = lines(
    "ind\tpop\tL1\tL1",
    "i1\tA\t10\t11",
    "i1\tA\t10\t12",
  );
  const g = parseStrafTsv(txt, 2);
  const w = getParseWarnings(g);
  assert.ok(w.some((i) => /Duplicate individual/i.test(i.message)));
});

// --- 8. Issue object structure --------------------------------------------

test("ParseError carries structured issues with line/column info", () => {
  const txt = lines("ind\tpop\tL1\tL1", "i1\tA\t10\tBAD");
  const issues = expectParseError(txt, 2);
  const err = issues.find((i) => /not numeric/i.test(i.message));
  assert.ok(err, "missing 'not numeric' error");
  assert.equal(err!.line, 2);
  assert.ok(err!.column !== undefined);
  assert.ok(err!.column! >= 1);
});

// --- 9. ParseError.message is a readable summary --------------------------

test("ParseError.message is a multi-line summary", () => {
  const txt = lines("foo\tbar\tL1", "i1\tA\t10");
  try {
    parseStrafTsv(txt, 1);
    assert.fail("expected ParseError");
  } catch (e) {
    assert.ok(e instanceof ParseError);
    assert.ok(e.message.includes("ERROR"));
    assert.ok(e.message.includes("ind"));
    assert.ok(e.message.includes("pop"));
  }
});
