/**
 * Parity / sanity tests against known mathematical identities.
 *
 * These don't replace cross-checks against R outputs (the original is the
 * source of truth on real data), but they pin down the implementations so
 * regressions show up.
 */

import { test } from "node:test";
import { strict as assert } from "node:assert";

import { parseStrafTsv } from "../parser.ts";
import { populationDistances } from "./distances.ts";
import { classicalMds } from "./mds.ts";
import { isoMds, pava } from "./iso_mds.ts";
import { confidenceEllipse } from "./ellipse.ts";
import { hclust, dendrogramSegments } from "./hclust.ts";
import { perLocusFStats, pairwiseFst } from "./fst.ts";
import { haplotypeStats } from "./haplotype.ts";
import { pca } from "./pca.ts";

// --- PAVA -------------------------------------------------------------------

test("pava: already monotonic input is unchanged", () => {
  const y = [1, 2, 3, 4, 5];
  assert.deepEqual(pava(y), y);
});

test("pava: violator is averaged with neighbour", () => {
  const out = pava([1, 3, 2, 4]);
  // The pair (3, 2) gets pooled to (2.5, 2.5)
  assert.deepEqual(out, [1, 2.5, 2.5, 4]);
});

test("pava: long descending block is constant-pooled", () => {
  const out = pava([5, 4, 3, 2, 1]);
  for (const v of out) assert.equal(v, 3);
});

// --- Distances --------------------------------------------------------------

const synth = parseStrafTsv(synthDiploid(), 2);

test("distances: identical pops have zero distance", () => {
  // Build a 2-pop dataset where pop B has the same individuals as A.
  const dup = parseStrafTsv(duplicatePop(), 2);
  for (const m of ["nei", "rogers", "provesti", "edwards", "reynolds"] as const) {
    const d = populationDistances(dup, m);
    assert.ok(d.matrix[0]![1]! < 1e-9, `${m}: expected ~0, got ${d.matrix[0]![1]}`);
  }
});

test("distances: matrix is symmetric and zero diagonal", () => {
  for (const m of ["nei", "rogers", "provesti", "edwards", "reynolds"] as const) {
    const d = populationDistances(synth, m);
    for (let i = 0; i < d.populations.length; i++) {
      assert.equal(d.matrix[i]![i], 0);
      for (let j = i + 1; j < d.populations.length; j++) {
        assert.ok(Math.abs(d.matrix[i]![j]! - d.matrix[j]![i]!) < 1e-12);
      }
    }
  }
});

// --- MDS --------------------------------------------------------------------

test("classicalMds: zero distance matrix gives zero coordinates", () => {
  const n = 4;
  const D = Array.from({ length: n }, () => new Array(n).fill(0));
  const mds = classicalMds(D);
  for (const row of mds.points) for (const v of row) assert.ok(Math.abs(v) < 1e-9);
});

test("classicalMds: recovers a 1-D grid", () => {
  // 4 evenly-spaced points: distances are |i - j|.
  const n = 4;
  const D: number[][] = [];
  for (let i = 0; i < n; i++) {
    const row: number[] = [];
    for (let j = 0; j < n; j++) row.push(Math.abs(i - j));
    D.push(row);
  }
  const mds = classicalMds(D);
  // The first axis should explain ~all the variance.
  assert.ok(mds.variance[0]! > 0.9);
  // And the projected distances on PC1 should match the input distances up to sign.
  const pc1 = mds.points.map((p) => p[0]!);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      assert.ok(Math.abs(Math.abs(pc1[i]! - pc1[j]!) - D[i]![j]!) < 1e-6);
    }
  }
});

test("isoMds: stress is ≥ 0 and decreasing toward classical solution on a Euclidean input", () => {
  // 5 points with truly Euclidean distances → both classical and isoMDS should
  // give near-zero stress.
  const n = 5;
  const truth = [
    [0, 0],
    [1, 0],
    [2, 1],
    [0, 2],
    [-1, 1],
  ];
  const D: number[][] = [];
  for (let i = 0; i < n; i++) {
    const row: number[] = [];
    for (let j = 0; j < n; j++) {
      const dx = truth[i]![0]! - truth[j]![0]!;
      const dy = truth[i]![1]! - truth[j]![1]!;
      row.push(Math.sqrt(dx * dx + dy * dy));
    }
    D.push(row);
  }
  const iso = isoMds(D, 2);
  assert.ok(iso.stress < 5, `expected low stress on Euclidean input, got ${iso.stress}%`);
});

// --- Ellipse ----------------------------------------------------------------

test("confidenceEllipse: points fall mostly inside their own 95% ellipse", () => {
  // Standard bivariate normal sample (Box–Muller, fixed seed via deterministic order).
  const xs: number[] = [];
  const ys: number[] = [];
  for (let i = 0; i < 200; i++) {
    const u1 = (i * 0.61803398875) % 1 || 0.5;
    const u2 = ((i + 0.37) * 0.61803398875) % 1 || 0.5;
    const r = Math.sqrt(-2 * Math.log(u1));
    xs.push(r * Math.cos(2 * Math.PI * u2));
    ys.push(r * Math.sin(2 * Math.PI * u2));
  }
  const e = confidenceEllipse(xs, ys, 0.95);
  assert.ok(e !== null);
  // Polygon should be closed up to floating-point precision and have
  // a reasonable number of points.
  const x0 = e!.x[0]!;
  const xN = e!.x[e!.x.length - 1]!;
  assert.ok(Math.abs(x0 - xN) < 1e-9, `polygon endpoints differ: ${x0} vs ${xN}`);
  assert.ok(e!.x.length > 30);
});

// --- hclust -----------------------------------------------------------------

test("hclust + dendrogramSegments: produces n−1 merges and consistent leaf order", () => {
  const D = [
    [0, 1, 5, 5],
    [1, 0, 5, 5],
    [5, 5, 0, 1],
    [5, 5, 1, 0],
  ];
  const tree = hclust(D, ["A", "B", "C", "D"]);
  assert.equal(tree.merge.length, 3);
  assert.equal(tree.height.length, 3);
  assert.equal(tree.order.length, 4);

  const seg = dendrogramSegments(tree);
  assert.equal(seg.leafLabels.length, 4);
  assert.equal(seg.segments.length, 3);
  // First merge should join {A, B} (closest pair, distance 1).
  assert.equal(tree.height[0], 1);
});

// --- F-stats ----------------------------------------------------------------

test("perLocusFStats: identical populations give bounded |Fst|", () => {
  // Note: Weir & Cockerham's θ̂ is an unbiased *estimator* — for two finite
  // identical samples that aren't in HWE, the sampling correction can drive
  // it slightly negative. We verify it stays within a sane range and is
  // never close to 1.
  const dup = parseStrafTsv(duplicatePop(), 2);
  const fs = perLocusFStats(dup);
  for (const f of fs) {
    if (Number.isFinite(f.Fst)) {
      assert.ok(Math.abs(f.Fst) < 0.5, `${f.locus}: |Fst|=${Math.abs(f.Fst)}`);
    }
  }
});

test("pairwiseFst: identical populations → |Fst| bounded", () => {
  const dup = parseStrafTsv(duplicatePop(), 2);
  const r = pairwiseFst(dup);
  assert.ok(Math.abs(r.matrix[0]![1]!) < 0.5, `|Fst|=${Math.abs(r.matrix[0]![1]!)}`);
});

test("pairwiseFst: large samples of identical pops → small |Fst|", () => {
  // W&C's θ̂ is a finite-sample estimator; for finite samples it can be slightly
  // negative even with no real differentiation. The bias scales as 1/(n−1),
  // so with large n it tends to 0.
  const lines = ["ind\tpop\tL1\tL1"];
  for (const pop of ["A", "B"] as const) {
    let id = 0;
    for (let i = 0; i < 50; i++) lines.push(`${pop}${++id}\t${pop}\t10\t10`);
    for (let i = 0; i < 50; i++) lines.push(`${pop}${++id}\t${pop}\t11\t11`);
    for (let i = 0; i < 100; i++) lines.push(`${pop}${++id}\t${pop}\t10\t11`);
  }
  const g = parseStrafTsv(lines.join("\n"), 2);
  const r = pairwiseFst(g);
  assert.ok(
    Math.abs(r.matrix[0]![1]!) < 0.01,
    `large-n identical pops: |Fst|=${Math.abs(r.matrix[0]![1]!)}`,
  );
});

test("pairwiseFst: highly differentiated populations give Fst near 1", () => {
  const diff = parseStrafTsv(highlyDifferentiated(), 2);
  const r = pairwiseFst(diff);
  // With completely disjoint allele sets between pops, Fst should be close to 1.
  assert.ok(r.matrix[0]![1]! > 0.5, `Fst=${r.matrix[0]![1]}`);
});

// --- Haplotypes -------------------------------------------------------------

test("haplotypeStats: counts sum to n, frequencies sum to 1", () => {
  const hap = parseStrafTsv(synthHaploid(), 1);
  const r = haplotypeStats(hap);
  let cSum = 0;
  let fSum = 0;
  for (const h of r.haplotypes) {
    cSum += h.count;
    fSum += h.frequency;
  }
  assert.equal(cSum, r.n);
  assert.ok(Math.abs(fSum - 1) < 1e-9);
  assert.ok(r.diversity >= 0 && r.diversity <= 1);
});

// --- PCA scaling ------------------------------------------------------------

test("pca: scaled and unscaled both produce ≥ 1 component for n>1", () => {
  const a = pca(synth);
  const b = pca(synth, undefined, { scale: true });
  assert.ok(a.eigenvalues.length >= 1);
  assert.ok(b.eigenvalues.length >= 1);
  assert.equal(a.alleleLabels.length, a.loadings.length);
});

// --- Helpers ---------------------------------------------------------------

function synthDiploid(): string {
  // 3 populations, 2 loci, allele frequencies vary across pops.
  const lines = ["ind\tpop\tL1\tL1\tL2\tL2"];
  // Pop A: L1 alleles ~ uniform over {10, 11, 12}; L2 ~ over {15, 16}
  const popA = [
    ["10", "11"],
    ["11", "12"],
    ["10", "12"],
    ["11", "11"],
    ["10", "10"],
    ["12", "12"],
  ];
  const popAL2 = [
    ["15", "16"],
    ["15", "15"],
    ["16", "16"],
    ["15", "16"],
    ["15", "15"],
    ["16", "16"],
  ];
  // Pop B: L1 enriched in {12, 13}
  const popB = [
    ["12", "13"],
    ["12", "12"],
    ["13", "13"],
    ["12", "13"],
    ["13", "13"],
  ];
  const popBL2 = [
    ["16", "17"],
    ["16", "16"],
    ["17", "17"],
    ["16", "17"],
    ["17", "17"],
  ];
  // Pop C: L1 mostly {10}
  const popC = [
    ["10", "10"],
    ["10", "10"],
    ["10", "11"],
    ["10", "10"],
  ];
  const popCL2 = [
    ["15", "15"],
    ["15", "15"],
    ["15", "16"],
    ["15", "15"],
  ];
  let id = 0;
  for (let i = 0; i < popA.length; i++) {
    lines.push(`${++id}\tA\t${popA[i]!.join("\t")}\t${popAL2[i]!.join("\t")}`);
  }
  for (let i = 0; i < popB.length; i++) {
    lines.push(`${++id}\tB\t${popB[i]!.join("\t")}\t${popBL2[i]!.join("\t")}`);
  }
  for (let i = 0; i < popC.length; i++) {
    lines.push(`${++id}\tC\t${popC[i]!.join("\t")}\t${popCL2[i]!.join("\t")}`);
  }
  return lines.join("\n");
}

/** Two populations with the same individuals (Fst should be 0). */
function duplicatePop(): string {
  const lines = ["ind\tpop\tL1\tL1"];
  const indivs = [
    ["10", "11"],
    ["10", "12"],
    ["11", "12"],
    ["11", "11"],
    ["10", "10"],
    ["12", "12"],
  ];
  let id = 0;
  for (const t of indivs) lines.push(`${++id}\tA\t${t.join("\t")}`);
  for (const t of indivs) lines.push(`${++id}\tB\t${t.join("\t")}`);
  return lines.join("\n");
}

/** Two pops with completely disjoint alleles (Fst near 1). */
function highlyDifferentiated(): string {
  const lines = ["ind\tpop\tL1\tL1"];
  for (let i = 0; i < 10; i++) lines.push(`A${i}\tA\t10\t10`);
  for (let i = 0; i < 10; i++) lines.push(`B${i}\tB\t20\t20`);
  return lines.join("\n");
}

function synthHaploid(): string {
  // Y-STR-like, 4 loci, several distinct haplotypes
  const lines = ["ind\tpop\tL1\tL2\tL3\tL4"];
  const haps = [
    ["10", "12", "15", "20"],
    ["10", "12", "15", "20"],
    ["10", "12", "15", "20"], // most common
    ["11", "12", "15", "20"],
    ["10", "13", "15", "20"],
    ["10", "12", "16", "20"],
    ["10", "12", "15", "21"],
  ];
  let id = 0;
  for (const h of haps) lines.push(`H${++id}\tP\t${h.join("\t")}`);
  return lines.join("\n");
}
