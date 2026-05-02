/**
 * Classical multidimensional scaling (a.k.a. principal coordinates analysis,
 * a.k.a. Torgerson scaling) on a square distance matrix.
 *
 * Algorithm:
 *   1. Square the distance matrix elementwise to get D².
 *   2. Double-center: B = -½ · J · D² · J  where J = I − (1/n) · 11ᵀ
 *      Equivalently B_ij = -½ (D²_ij − rowMean_i − colMean_j + grandMean)
 *   3. Eigendecompose B; coordinates = eigenvectors · sqrt(max(0, λ))
 *
 * Eigendecomposition uses the same Jacobi routine as PCA. Negative eigenvalues
 * (which can arise for non-Euclidean distances) are clamped to zero.
 */

export interface MdsResult {
  /** points[i][k] = coordinate of object i on axis k */
  points: number[][];
  eigenvalues: number[];
  /** Proportion of variance explained per axis (negative eigenvalues excluded). */
  variance: number[];
}

export function classicalMds(distance: number[][]): MdsResult {
  const n = distance.length;
  if (n < 2) throw new Error("MDS requires at least 2 objects.");

  // Squared distances.
  const D2: number[][] = distance.map((row) => row.map((d) => d * d));

  // Row means, column means, grand mean.
  const rowMean = new Array(n).fill(0) as number[];
  const colMean = new Array(n).fill(0) as number[];
  let grandMean = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const v = D2[i]![j]!;
      rowMean[i]! += v;
      colMean[j]! += v;
      grandMean += v;
    }
  }
  for (let i = 0; i < n; i++) {
    rowMean[i] = rowMean[i]! / n;
    colMean[i] = colMean[i]! / n;
  }
  grandMean /= n * n;

  const B: number[][] = Array.from({ length: n }, () => new Array(n).fill(0) as number[]);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      B[i]![j] = -0.5 * (D2[i]![j]! - rowMean[i]! - colMean[j]! + grandMean);
    }
  }

  const { values, vectors } = jacobiEigenSymmetric(B);

  // Coordinates: V * diag(sqrt(max(0, λ))). Vectors come back as columns.
  const points: number[][] = Array.from({ length: n }, () => new Array(n).fill(0) as number[]);
  for (let i = 0; i < n; i++) {
    for (let k = 0; k < n; k++) {
      points[i]![k] = vectors[i]![k]! * Math.sqrt(Math.max(0, values[k]!));
    }
  }

  const positiveSum = values.reduce((s, v) => s + Math.max(0, v), 0);
  const variance = values.map((v) =>
    positiveSum > 0 ? Math.max(0, v) / positiveSum : 0,
  );

  return { points, eigenvalues: values, variance };
}

/**
 * Jacobi eigendecomposition for a real symmetric matrix.
 * Returns eigenvalues sorted descending and matching eigenvectors as columns.
 *
 * (Duplicated from pca.ts to keep this module self-contained — both copies
 * are small and behavior is identical.)
 */
function jacobiEigenSymmetric(A: number[][]): { values: number[]; vectors: number[][] } {
  const n = A.length;
  const a: number[][] = Array.from({ length: n }, (_, i) => A[i]!.slice());
  const v: number[][] = Array.from({ length: n }, (_, i) => {
    const row = new Array(n).fill(0) as number[];
    row[i] = 1;
    return row;
  });

  for (let sweep = 0; sweep < 100; sweep++) {
    let off = 0;
    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) off += Math.abs(a[p]![q]!);
    }
    if (off < 1e-12 * n * n) break;

    for (let p = 0; p < n - 1; p++) {
      for (let q = p + 1; q < n; q++) {
        const apq = a[p]![q]!;
        if (Math.abs(apq) < 1e-14) continue;
        const app = a[p]![p]!;
        const aqq = a[q]![q]!;
        const theta = (aqq - app) / (2 * apq);
        let t: number;
        if (Math.abs(theta) > 1e15) t = 1 / (2 * theta);
        else t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1);
        const s = t * c;

        a[p]![p] = app - t * apq;
        a[q]![q] = aqq + t * apq;
        a[p]![q] = 0;
        a[q]![p] = 0;
        for (let r = 0; r < n; r++) {
          if (r !== p && r !== q) {
            const arp = a[r]![p]!;
            const arq = a[r]![q]!;
            a[r]![p] = c * arp - s * arq;
            a[p]![r] = a[r]![p]!;
            a[r]![q] = s * arp + c * arq;
            a[q]![r] = a[r]![q]!;
          }
          const vrp = v[r]![p]!;
          const vrq = v[r]![q]!;
          v[r]![p] = c * vrp - s * vrq;
          v[r]![q] = s * vrp + c * vrq;
        }
      }
    }
  }

  const values = new Array(n).fill(0).map((_, i) => a[i]![i]!);
  const order = values.map((_, i) => i).sort((i, j) => values[j]! - values[i]!);
  const sortedValues = order.map((i) => values[i]!);
  const sortedVectors: number[][] = Array.from({ length: n }, () => new Array(n).fill(0) as number[]);
  for (let r = 0; r < n; r++) {
    for (let c = 0; c < n; c++) sortedVectors[r]![c] = v[r]![order[c]!]!;
  }
  return { values: sortedValues, vectors: sortedVectors };
}
