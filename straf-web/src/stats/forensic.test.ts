import { test } from "node:test";
import { strict as assert } from "node:assert";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, resolve } from "node:path";

import { parseStrafTsv } from "../parser.ts";
import { locusIndices } from "./forensic.ts";
import { hweChiSquare } from "./hwe.ts";
import { alleleFrequencies } from "./freq.ts";

const here = dirname(fileURLToPath(import.meta.url));
const examplePath = resolve(here, "../../public/exampleSTRAFdiplo.txt");
const text = readFileSync(examplePath, "utf-8");
const genos = parseStrafTsv(text, 2);

test("parser: shape", () => {
  assert.equal(genos.ploidy, 2);
  assert.ok(genos.individuals.length > 0);
  assert.ok(genos.loci.length > 0);
  assert.equal(genos.alleles.length, genos.loci.length);
  for (const col of genos.alleles) assert.equal(col.length, genos.individuals.length);
});

test("alleleFrequencies: per-locus columns sum to 1 (excluding missing)", () => {
  const tbl = alleleFrequencies(genos);
  for (let l = 0; l < tbl.loci.length; l++) {
    let sum = 0;
    let any = false;
    for (let a = 0; a < tbl.alleles.length; a++) {
      const v = tbl.values[a]![l];
      if (v !== null && v !== undefined) {
        sum += v;
        any = true;
      }
    }
    if (any) assert.ok(Math.abs(sum - 1) < 1e-9, `locus ${tbl.loci[l]}: freq sum = ${sum}`);
  }
});

test("forensic: GD and PIC are in (0,1) for polymorphic STR loci", () => {
  const stats = locusIndices(genos);
  for (const s of stats) {
    assert.ok(s.GD > 0 && s.GD < 1, `${s.locus} GD=${s.GD}`);
    assert.ok(s.PIC > 0 && s.PIC < 1, `${s.locus} PIC=${s.PIC}`);
    assert.ok(s.PIC <= s.GD + 1e-9, `${s.locus} PIC > GD`);
    assert.ok(s.PD >= 0 && s.PD <= 1);
    assert.ok(s.PM >= 0 && s.PM <= 1);
    assert.ok(Math.abs(s.PM + s.PD - 1) < 1e-9);
    if (s.Hobs !== undefined) {
      assert.ok(s.Hobs >= 0 && s.Hobs <= 1);
      assert.ok(s.PE !== undefined);
      assert.ok(s.TPI !== undefined);
    }
  }
});

test("forensic: D3S1358 has expected number of alleles", () => {
  const stats = locusIndices(genos);
  const d3 = stats.find((s) => s.locus === "D3S1358");
  assert.ok(d3, "D3S1358 not found");
  // Spot check: this dataset is small but D3S1358 is highly polymorphic
  assert.ok(d3!.Nall >= 4, `D3S1358 Nall=${d3!.Nall}`);
});

test("hwe: returns one result per locus, p-values in [0, 1] or NaN", () => {
  const res = hweChiSquare(genos);
  assert.equal(res.length, genos.loci.length);
  for (const r of res) {
    if (Number.isFinite(r.pValue)) {
      assert.ok(r.pValue >= 0 && r.pValue <= 1, `${r.locus} p=${r.pValue}`);
      assert.ok(r.chiSq >= 0);
      assert.ok(r.df > 0);
    }
  }
});

test("hwe: chi-square monotonicity sanity (large chi2 => small p)", () => {
  // Construct a tiny synthetic dataset that strongly violates HWE
  // (all heterozygotes, no homozygotes — Wahlund-like deviation).
  const synthHeader = "ind\tpop\tL1\tL1";
  const rows: string[] = [synthHeader];
  for (let i = 0; i < 100; i++) rows.push(`I${i}\tP\t10\t12`);
  const synth = parseStrafTsv(rows.join("\n"), 2);
  const r = hweChiSquare(synth)[0]!;
  assert.ok(r.pValue < 0.05, `expected p<0.05 for all-heterozygous synthetic, got ${r.pValue}`);
});
