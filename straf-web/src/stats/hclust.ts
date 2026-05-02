/**
 * Hierarchical agglomerative clustering on a square distance matrix.
 *
 * Default linkage is UPGMA (average), matching `stats::hclust(dst, method="average")`
 * which is the standard choice for genetic distance trees in popgen — and what
 * `ape::as.phylo(hclust(...))` produces in the original STRAF.
 *
 * Output is the same shape as R's hclust: `merge` is an (n−1) × 2 integer
 * matrix (negative IDs are leaves, positive are previous merges, 1-indexed),
 * `height` is the merge height of each step, and `order` is a leaf permutation
 * for non-crossing tree drawing.
 */

export type Linkage = "single" | "complete" | "average";

export interface HclustResult {
  /** merge[k] = [a, b], 1-indexed; negatives = leaves, positives = previous merges. */
  merge: [number, number][];
  /** Merge heights, length n−1. */
  height: number[];
  /** Leaf order (0-indexed) for non-crossing tree layout. */
  order: number[];
  /** Original labels (passed through, unchanged). */
  labels: string[];
}

export function hclust(
  distance: number[][],
  labels: string[],
  linkage: Linkage = "average",
): HclustResult {
  const n = distance.length;
  if (n < 2) throw new Error("hclust requires at least 2 objects.");

  // Distance matrix copy — we'll update it in place as clusters merge.
  const D: number[][] = distance.map((row) => row.slice());

  // Track active clusters by id. Negative ids are leaves (-1..-n);
  // positive ids are previous merges (1..n-1). Use an ordered list of (id, size).
  type Active = { id: number; size: number; index: number /* row/col in D */ };
  const active: Active[] = [];
  for (let i = 0; i < n; i++) active.push({ id: -(i + 1), size: 1, index: i });

  const merge: [number, number][] = [];
  const height: number[] = [];
  // For drawing order we keep track of leaf sequences inside each cluster.
  const leafSeq: Map<number, number[]> = new Map();
  for (let i = 0; i < n; i++) leafSeq.set(-(i + 1), [i]);

  for (let step = 0; step < n - 1; step++) {
    // Find the closest pair among active clusters.
    let best = Infinity;
    let bi = -1;
    let bj = -1;
    for (let i = 0; i < active.length - 1; i++) {
      for (let j = i + 1; j < active.length; j++) {
        const d = D[active[i]!.index]![active[j]!.index]!;
        if (d < best) {
          best = d;
          bi = i;
          bj = j;
        }
      }
    }
    if (bi < 0 || bj < 0) throw new Error("hclust: no pair found");

    const ai = active[bi]!;
    const aj = active[bj]!;
    const newId = step + 1;
    merge.push([ai.id, aj.id]);
    height.push(best);
    // Merge leaf sequences (preserves insertion order → non-crossing layout).
    const seq = [...(leafSeq.get(ai.id) ?? []), ...(leafSeq.get(aj.id) ?? [])];
    leafSeq.set(newId, seq);

    // Update distance from the new cluster to every other active cluster.
    // We reuse one of the freed rows/cols (ai.index) for the new cluster.
    const newIdx = ai.index;
    for (let k = 0; k < active.length; k++) {
      if (k === bi || k === bj) continue;
      const idxK = active[k]!.index;
      const dik = D[ai.index]![idxK]!;
      const djk = D[aj.index]![idxK]!;
      let dnk: number;
      switch (linkage) {
        case "single":
          dnk = Math.min(dik, djk);
          break;
        case "complete":
          dnk = Math.max(dik, djk);
          break;
        case "average": {
          const wi = ai.size;
          const wj = aj.size;
          dnk = (wi * dik + wj * djk) / (wi + wj);
          break;
        }
      }
      D[newIdx]![idxK] = dnk;
      D[idxK]![newIdx] = dnk;
    }
    D[newIdx]![newIdx] = 0;

    // Remove aj from active, replace ai with the merged cluster.
    const merged: Active = { id: newId, size: ai.size + aj.size, index: newIdx };
    active.splice(bj, 1);
    active.splice(bi, 1, merged);
  }

  const finalId = active[0]!.id;
  const order = leafSeq.get(finalId) ?? [];
  return { merge, height, order, labels };
}

/**
 * Convert hclust output into line segments suitable for plotting a horizontal
 * dendrogram (root on top, leaves on bottom).
 *
 * Returns:
 *   - leafX[i]    = x position of the i-th leaf (in the original label order)
 *   - segments    = list of polylines (each is an {x, y} pair list) describing
 *                   the U-shaped joins; render them as a plotly scatter trace
 *                   with mode="lines"
 *   - leafLabels  = labels in their plotting order (left to right)
 */
export function dendrogramSegments(h: HclustResult): {
  leafX: number[];
  leafLabels: string[];
  segments: { x: number[]; y: number[] }[];
  height: number;
} {
  const n = h.labels.length;
  const orderPos = new Array<number>(n).fill(0);
  for (let i = 0; i < h.order.length; i++) orderPos[h.order[i]!] = i;

  // For each cluster id, store its anchor (x, y).
  // Leaves sit at y = 0; merges sit at y = merge height.
  const anchor = new Map<number, { x: number; y: number }>();
  for (let i = 0; i < n; i++) {
    anchor.set(-(i + 1), { x: orderPos[i]!, y: 0 });
  }

  const segments: { x: number[]; y: number[] }[] = [];
  for (let s = 0; s < h.merge.length; s++) {
    const [a, b] = h.merge[s]!;
    const A = anchor.get(a)!;
    const B = anchor.get(b)!;
    const y = h.height[s]!;
    const xMid = (A.x + B.x) / 2;

    // U shape: down from (A.x, y) → (A.x, A.y), across (A.x, y) → (B.x, y),
    // down (B.x, y) → (B.x, B.y).
    segments.push({ x: [A.x, A.x, B.x, B.x], y: [A.y, y, y, B.y] });
    anchor.set(s + 1, { x: xMid, y });
  }

  return {
    leafX: h.order.map((_, i) => i),
    leafLabels: h.order.map((leaf) => h.labels[leaf]!),
    segments,
    height: h.height[h.height.length - 1] ?? 0,
  };
}
