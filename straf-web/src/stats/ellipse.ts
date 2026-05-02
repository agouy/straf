/**
 * 95% (or arbitrary level) confidence ellipse for 2D point cloud, assuming
 * bivariate normality — same convention as ggplot2's `stat_ellipse(type = "norm")`.
 *
 * Returns a polygon (closed) approximation of the ellipse boundary.
 *
 * Math: with sample mean μ and sample covariance Σ, the level-α ellipse is
 *   { x : (x − μ)ᵀ Σ⁻¹ (x − μ) = q }
 * where q = χ²(2 dof, α). For α=0.95 → q = 5.991.
 *
 * We diagonalize Σ = R diag(λ) Rᵀ, then the ellipse is parameterized as
 *   x(t) = μ + sqrt(q) · R · diag(sqrt(λ)) · (cos t, sin t)ᵀ
 */
export function confidenceEllipse(
  xs: number[],
  ys: number[],
  level = 0.95,
  nPoints = 64,
): { x: number[]; y: number[] } | null {
  const n = xs.length;
  if (n < 3) return null;

  let mx = 0;
  let my = 0;
  for (let i = 0; i < n; i++) {
    mx += xs[i]!;
    my += ys[i]!;
  }
  mx /= n;
  my /= n;

  let sxx = 0;
  let syy = 0;
  let sxy = 0;
  for (let i = 0; i < n; i++) {
    const dx = xs[i]! - mx;
    const dy = ys[i]! - my;
    sxx += dx * dx;
    syy += dy * dy;
    sxy += dx * dy;
  }
  const denom = n - 1;
  if (denom <= 0) return null;
  sxx /= denom;
  syy /= denom;
  sxy /= denom;

  // 2x2 eigendecomposition (closed-form).
  const tr = sxx + syy;
  const det = sxx * syy - sxy * sxy;
  const disc = Math.sqrt(Math.max(0, tr * tr / 4 - det));
  const lam1 = tr / 2 + disc;
  const lam2 = tr / 2 - disc;
  if (lam1 <= 0 || lam2 < 0) return null;

  // Eigenvectors. If sxy is non-zero, eigenvector for lam1 is (lam1 - syy, sxy)
  // (proportional). Otherwise the matrix is diagonal.
  let v1x: number;
  let v1y: number;
  if (Math.abs(sxy) > 1e-12) {
    v1x = lam1 - syy;
    v1y = sxy;
  } else if (sxx >= syy) {
    v1x = 1;
    v1y = 0;
  } else {
    v1x = 0;
    v1y = 1;
  }
  const norm1 = Math.hypot(v1x, v1y);
  v1x /= norm1;
  v1y /= norm1;
  const v2x = -v1y;
  const v2y = v1x;

  const q = chiSquareQuantileDof2(level);
  const a = Math.sqrt(q * lam1);
  const b = Math.sqrt(q * lam2);

  const xOut: number[] = [];
  const yOut: number[] = [];
  for (let i = 0; i <= nPoints; i++) {
    const t = (2 * Math.PI * i) / nPoints;
    const c = Math.cos(t);
    const s = Math.sin(t);
    const x = mx + a * c * v1x + b * s * v2x;
    const y = my + a * c * v1y + b * s * v2y;
    xOut.push(x);
    yOut.push(y);
  }
  return { x: xOut, y: yOut };
}

/**
 * Quantile function of chi-square distribution with 2 degrees of freedom:
 * F(x) = 1 − exp(−x/2), so F⁻¹(p) = −2 ln(1 − p).
 * No iteration needed.
 */
function chiSquareQuantileDof2(p: number): number {
  if (p <= 0) return 0;
  if (p >= 1) return Infinity;
  return -2 * Math.log(1 - p);
}
