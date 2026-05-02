/**
 * Kruskal's non-metric multidimensional scaling (isoMDS).
 *
 * Approach: SMACOF-style majorization on Kruskal's stress, with monotonic
 * regression of fitted distances against the input dissimilarities at every
 * iteration. This is the same shape of algorithm as `MASS::isoMDS`, though
 * MASS uses a numerical optimizer; SMACOF is more numerically stable and
 * trivially convergent.
 *
 *   stress = sqrt( Σ (d_ij − δ̂_ij)² / Σ d_ij² )
 *
 * where d_ij are Euclidean distances in the embedding and δ̂_ij are the
 * monotonic regression of d_ij against the input dissimilarities — the
 * "disparities". PAVA (pool-adjacent-violators) computes the monotonic
 * regression in O(n log n).
 *
 * Initial configuration is taken from classical MDS (cmdscale) to give the
 * iteration a sane starting point — same convention as MASS::isoMDS.
 */

import { classicalMds } from "./mds.ts";

export interface IsoMdsResult {
  /** points[i][k] = coordinate of object i on axis k */
  points: number[][];
  /** Final Kruskal stress (×100 to match MASS::isoMDS reporting). */
  stress: number;
  iterations: number;
  converged: boolean;
}

export function isoMds(
  distance: number[][],
  k = 2,
  opts: { maxIter?: number; tol?: number } = {},
): IsoMdsResult {
  const n = distance.length;
  if (n < 3) throw new Error("isoMDS requires at least 3 objects.");
  if (k < 1 || k > n - 1) throw new Error(`k=${k} out of range for n=${n}.`);

  const maxIter = opts.maxIter ?? 200;
  const tol = opts.tol ?? 1e-5;

  // Initialize from classical MDS, take first k columns.
  const cmds = classicalMds(distance);
  let X: number[][] = cmds.points.map((row) => row.slice(0, k));
  // Pad if k > available components.
  if (X[0]!.length < k) {
    for (const row of X) while (row.length < k) row.push(0);
  }

  // Pre-compute the ordering of pairs (i, j) by input dissimilarity ascending,
  // for monotonic regression. Ties are pooled in PAVA, so any stable order is fine.
  const pairs: { i: number; j: number; delta: number }[] = [];
  for (let i = 0; i < n - 1; i++) {
    for (let j = i + 1; j < n; j++) {
      pairs.push({ i, j, delta: distance[i]![j]! });
    }
  }
  pairs.sort((a, b) => a.delta - b.delta);

  let prevStress = Infinity;
  let stress = NaN;
  let iter = 0;
  let converged = false;

  for (iter = 0; iter < maxIter; iter++) {
    // 1. Compute fitted Euclidean distances d_ij in current configuration.
    const d: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let i = 0; i < n - 1; i++) {
      for (let j = i + 1; j < n; j++) {
        let s = 0;
        for (let a = 0; a < k; a++) {
          const diff = X[i]![a]! - X[j]![a]!;
          s += diff * diff;
        }
        const v = Math.sqrt(s);
        d[i]![j] = v;
        d[j]![i] = v;
      }
    }

    // 2. Monotonic regression (PAVA) of d_ij in the order of input δ_ij,
    //    yielding disparities δ̂_ij that are weakly increasing in δ.
    const dByPair = pairs.map((p) => d[p.i]![p.j]!);
    const disp = pava(dByPair);

    // Map back to a matrix.
    const Dhat: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let p = 0; p < pairs.length; p++) {
      const { i, j } = pairs[p]!;
      Dhat[i]![j] = disp[p]!;
      Dhat[j]![i] = disp[p]!;
    }

    // 3. Stress.
    let num = 0;
    let den = 0;
    for (let i = 0; i < n - 1; i++) {
      for (let j = i + 1; j < n; j++) {
        const dij = d[i]![j]!;
        const dhat = Dhat[i]![j]!;
        num += (dij - dhat) * (dij - dhat);
        den += dij * dij;
      }
    }
    stress = den > 0 ? Math.sqrt(num / den) : 0;

    if (Math.abs(prevStress - stress) < tol) {
      converged = true;
      break;
    }
    prevStress = stress;

    // 4. SMACOF-style update: X_new = (1/n) · B(X) · X
    //    B(X)_{ij} = -δ̂_ij / d_ij  for i ≠ j, d_ij > 0
    //    B(X)_{ii} = -Σ_{j ≠ i} B(X)_{ij}
    const B: number[][] = Array.from({ length: n }, () => new Array(n).fill(0));
    for (let i = 0; i < n; i++) {
      let rowSum = 0;
      for (let j = 0; j < n; j++) {
        if (i === j) continue;
        const dij = d[i]![j]!;
        const v = dij > 1e-12 ? -Dhat[i]![j]! / dij : 0;
        B[i]![j] = v;
        rowSum += v;
      }
      B[i]![i] = -rowSum;
    }

    const Xnew: number[][] = Array.from({ length: n }, () => new Array(k).fill(0));
    for (let i = 0; i < n; i++) {
      for (let a = 0; a < k; a++) {
        let s = 0;
        for (let j = 0; j < n; j++) s += B[i]![j]! * X[j]![a]!;
        Xnew[i]![a] = s / n;
      }
    }
    X = Xnew;
  }

  return {
    points: X,
    stress: stress * 100, // %, matching MASS::isoMDS convention
    iterations: iter,
    converged,
  };
}

/**
 * Pool-Adjacent-Violators Algorithm — monotonic non-decreasing regression.
 * Given y_1, ..., y_n in some natural order (here: order of input dissimilarities),
 * returns ŷ_1, ..., ŷ_n minimizing Σ (y_i − ŷ_i)² subject to ŷ_1 ≤ ŷ_2 ≤ ... ≤ ŷ_n.
 *
 * Standard implementation: walk left to right, merging blocks whose mean
 * violates the monotonicity constraint. O(n) amortized.
 */
export function pava(y: number[]): number[] {
  const n = y.length;
  if (n === 0) return [];

  // Each block: [startIndex, length, mean]
  const stack: { start: number; len: number; mean: number }[] = [];
  for (let i = 0; i < n; i++) {
    let cur = { start: i, len: 1, mean: y[i]! };
    while (stack.length > 0 && stack[stack.length - 1]!.mean > cur.mean) {
      const top = stack.pop()!;
      const newLen = top.len + cur.len;
      const newMean = (top.mean * top.len + cur.mean * cur.len) / newLen;
      cur = { start: top.start, len: newLen, mean: newMean };
    }
    stack.push(cur);
  }

  const out = new Array(n).fill(0) as number[];
  for (const block of stack) {
    for (let k = 0; k < block.len; k++) out[block.start + k] = block.mean;
  }
  return out;
}
