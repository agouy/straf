import type { Genotypes } from "../parser.ts";

/**
 * PCA of the individual-by-allele dosage matrix.
 *
 * The matrix X has one row per individual and one column per (locus, allele)
 * pair. Entries are dosage counts (0, 1, or 2 for diploid; 0 or 1 for haploid).
 * Missing alleles at a locus are mean-imputed using that locus's column means
 * before centering — this matches the convention used by `adegenet::scaleGen`
 * with `NA.method = "mean"`.
 *
 * Columns are centered (mean = 0) but not scaled; this is the most common
 * default in popgen PCA and matches the original STRAF behavior.
 *
 * Eigendecomposition uses Jacobi rotation on the smaller of X X^T or X^T X.
 * For typical STR datasets (n_ind < ~500, n_alleles < ~300) this runs in
 * milliseconds.
 */
export interface PcaResult {
  /** scores[indIdx][componentIdx]; PC1 first, ordered by descending eigenvalue. */
  scores: number[][];
  /** eigenvalues (variance explained, in descending order). */
  eigenvalues: number[];
  /** Proportion of total variance explained per component. */
  variance: number[];
  /** Cumulative variance proportion. */
  cumulative: number[];
  /** Loadings (rotation matrix): rows = original allele columns, cols = PCs. */
  loadings: number[][];
  /** Allele-column labels in `<locus>.<allele>` form for the loadings rows. */
  alleleLabels: string[];
  individuals: string[];
  populations: string[];
}

export interface PcaOptions {
  /** Scale columns to unit variance after centering (matches prcomp(scale=TRUE)). Default: false. */
  scale?: boolean;
}

export function pca(
  genos: Genotypes,
  indMask?: boolean[],
  options: PcaOptions = {},
): PcaResult {
  const mask = indMask ?? genos.individuals.map(() => true);
  const indIdx: number[] = [];
  for (let i = 0; i < mask.length; i++) if (mask[i]) indIdx.push(i);
  const nInd = indIdx.length;
  if (nInd < 2) throw new Error("PCA requires at least 2 individuals.");

  // Build column index for every (locus, allele).
  const cols: { locus: number; allele: string }[] = [];
  for (let l = 0; l < genos.loci.length; l++) {
    const seen = new Set<string>();
    for (const t of genos.alleles[l]!) for (const a of t) if (a !== null) seen.add(a);
    const sorted = Array.from(seen).sort((a, b) => {
      const na = Number(a);
      const nb = Number(b);
      if (Number.isFinite(na) && Number.isFinite(nb)) return na - nb;
      return a.localeCompare(b);
    });
    for (const a of sorted) cols.push({ locus: l, allele: a });
  }
  const p = cols.length;

  // Build raw dosage matrix and a "missing" mask.
  const X: number[][] = Array.from({ length: nInd }, () => new Array(p).fill(0) as number[]);
  const missingByLocusInd: boolean[][] = Array.from({ length: nInd }, () =>
    new Array(genos.loci.length).fill(false) as boolean[],
  );

  for (let r = 0; r < nInd; r++) {
    const i = indIdx[r]!;
    for (let l = 0; l < genos.loci.length; l++) {
      const tuple = genos.alleles[l]![i]!;
      if (tuple.some((a) => a === null)) {
        missingByLocusInd[r]![l] = true;
      }
    }
  }

  // Fill dosages for non-missing genotypes.
  let cByLocusStart: number[] = new Array(genos.loci.length).fill(0);
  {
    let s = 0;
    for (let l = 0; l < genos.loci.length; l++) {
      cByLocusStart[l] = s;
      // count alleles for this locus
      let n = 0;
      for (const c of cols) if (c.locus === l) n++;
      s += n;
    }
  }
  const colByLocusAllele: Map<string, number> = new Map();
  cols.forEach((c, k) => colByLocusAllele.set(`${c.locus}|${c.allele}`, k));

  for (let r = 0; r < nInd; r++) {
    const i = indIdx[r]!;
    for (let l = 0; l < genos.loci.length; l++) {
      if (missingByLocusInd[r]![l]) continue;
      for (const a of genos.alleles[l]![i]!) {
        if (a === null) continue;
        const k = colByLocusAllele.get(`${l}|${a}`)!;
        X[r]![k]!++;
      }
    }
  }

  // Per-column means computed from non-missing rows only (locus-wise mask).
  const colMean = new Array(p).fill(0) as number[];
  const colDenom = new Array(p).fill(0) as number[];
  for (let r = 0; r < nInd; r++) {
    for (let k = 0; k < p; k++) {
      const l = cols[k]!.locus;
      if (missingByLocusInd[r]![l]) continue;
      colMean[k]! += X[r]![k]!;
      colDenom[k]!++;
    }
  }
  for (let k = 0; k < p; k++) {
    colMean[k] = colDenom[k]! > 0 ? colMean[k]! / colDenom[k]! : 0;
  }

  // Mean-impute missing rows then center.
  for (let r = 0; r < nInd; r++) {
    for (let k = 0; k < p; k++) {
      const l = cols[k]!.locus;
      if (missingByLocusInd[r]![l]) X[r]![k] = colMean[k]!;
      X[r]![k] = X[r]![k]! - colMean[k]!;
    }
  }

  // Optional scaling: divide each column by its standard deviation
  // (matching prcomp(scale=TRUE)). Columns with zero variance are left at 0.
  if (options.scale) {
    const colSd = new Array(p).fill(0) as number[];
    const denom = Math.max(1, nInd - 1);
    for (let k = 0; k < p; k++) {
      let s = 0;
      for (let r = 0; r < nInd; r++) s += X[r]![k]! * X[r]![k]!;
      colSd[k] = Math.sqrt(s / denom);
    }
    for (let r = 0; r < nInd; r++) {
      for (let k = 0; k < p; k++) {
        const sd = colSd[k]!;
        X[r]![k] = sd > 1e-12 ? X[r]![k]! / sd : 0;
      }
    }
  }

  // Eigendecompose. We always need loadings (eigenvectors of X^T X) for the
  // download — so compute them via the p × p path when feasible, or recover
  // them from the n × n decomposition when n < p.
  let eigenvalues: number[];
  let scores: number[][];
  let loadings: number[][]; // p × min(n, p)

  if (nInd <= p) {
    // n × n inner product — efficient when there are fewer individuals than allele columns
    const G = matMulNT(X, X); // X * X^T  : n × n
    const { values, vectors } = jacobiEigen(G);
    eigenvalues = values;
    scores = Array.from({ length: nInd }, () => new Array(values.length).fill(0));
    for (let r = 0; r < nInd; r++) {
      for (let c = 0; c < values.length; c++) {
        scores[r]![c] = vectors[r]![c]! * Math.sqrt(Math.max(0, values[c]!));
      }
    }
    // Recover loadings: V = X^T U / sqrt(λ)  (where U has columns of `vectors`)
    loadings = Array.from({ length: p }, () => new Array(values.length).fill(0));
    for (let c = 0; c < values.length; c++) {
      const sLam = Math.sqrt(Math.max(0, values[c]!));
      if (sLam < 1e-12) continue;
      for (let kk = 0; kk < p; kk++) {
        let s = 0;
        for (let r = 0; r < nInd; r++) s += X[r]![kk]! * vectors[r]![c]!;
        loadings[kk]![c] = s / sLam;
      }
    }
  } else {
    // p × p covariance form
    const C = matMulTN(X, X); // X^T * X : p × p
    const { values, vectors } = jacobiEigen(C);
    eigenvalues = values;
    scores = matMul(X, vectors);
    loadings = vectors;
  }

  // Variance explained: scale eigenvalues by 1/(n-1) for proportions, but the
  // *proportions* are unchanged.
  const total = eigenvalues.reduce((s, v) => s + Math.max(0, v), 0);
  const variance = eigenvalues.map((v) => (total > 0 ? Math.max(0, v) / total : 0));
  const cumulative: number[] = [];
  let acc = 0;
  for (const v of variance) {
    acc += v;
    cumulative.push(acc);
  }

  const alleleLabels = cols.map((c) => `${genos.loci[c.locus]}.${c.allele}`);

  return {
    scores,
    eigenvalues,
    variance,
    cumulative,
    loadings,
    alleleLabels,
    individuals: indIdx.map((i) => genos.individuals[i]!),
    populations: indIdx.map((i) => genos.populations[i]!),
  };
}

// --------- linear algebra helpers (small, allocation-light) ---------

function matMul(A: number[][], B: number[][]): number[][] {
  const n = A.length;
  const m = A[0]!.length;
  const p = B[0]!.length;
  const out: number[][] = Array.from({ length: n }, () => new Array(p).fill(0) as number[]);
  for (let i = 0; i < n; i++) {
    const Ai = A[i]!;
    const Oi = out[i]!;
    for (let k = 0; k < m; k++) {
      const a = Ai[k]!;
      if (a === 0) continue;
      const Bk = B[k]!;
      for (let j = 0; j < p; j++) Oi[j]! += a * Bk[j]!;
    }
  }
  return out;
}

/** Returns A * A^T (n × n, given A is n × m). */
function matMulNT(A: number[][], _A2: number[][]): number[][] {
  const n = A.length;
  const m = A[0]!.length;
  const out: number[][] = Array.from({ length: n }, () => new Array(n).fill(0) as number[]);
  for (let i = 0; i < n; i++) {
    const Ai = A[i]!;
    for (let j = i; j < n; j++) {
      const Aj = A[j]!;
      let s = 0;
      for (let k = 0; k < m; k++) s += Ai[k]! * Aj[k]!;
      out[i]![j] = s;
      out[j]![i] = s;
    }
  }
  return out;
}

/** Returns A^T * A (m × m). */
function matMulTN(A: number[][], _A2: number[][]): number[][] {
  const n = A.length;
  const m = A[0]!.length;
  const out: number[][] = Array.from({ length: m }, () => new Array(m).fill(0) as number[]);
  for (let i = 0; i < m; i++) {
    for (let j = i; j < m; j++) {
      let s = 0;
      for (let r = 0; r < n; r++) s += A[r]![i]! * A[r]![j]!;
      out[i]![j] = s;
      out[j]![i] = s;
    }
  }
  return out;
}

/**
 * Jacobi eigendecomposition for a real symmetric matrix.
 * Returns eigenvalues sorted descending and the matching eigenvectors as columns.
 */
function jacobiEigen(A: number[][]): { values: number[]; vectors: number[][] } {
  const n = A.length;
  // Copy A
  const a: number[][] = Array.from({ length: n }, (_, i) => A[i]!.slice());
  const v: number[][] = Array.from({ length: n }, (_, i) => {
    const row = new Array(n).fill(0) as number[];
    row[i] = 1;
    return row;
  });

  const maxSweeps = 100;
  for (let sweep = 0; sweep < maxSweeps; sweep++) {
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
        if (Math.abs(theta) > 1e15) {
          t = 1 / (2 * theta);
        } else {
          t = Math.sign(theta) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          if (theta === 0) t = 1;
        }
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
