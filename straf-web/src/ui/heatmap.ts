/**
 * Tiny SVG heatmap renderer.
 *
 * Used for plots where Observable Plot's marks aren't a good fit (LD p-value
 * matrix, haplotype pairwise differences). Hand-rolled because we want a
 * value rendered inside every cell — same behavior as the original R version.
 *
 * Hooks into the shared toolbar (`plotControls.ts`) so the heatmap also gets
 * a gear popover and SVG/PNG download buttons.
 */

import { attachToolbar, getFigureSlot } from "./plotControls.ts";
import type { ControlField } from "./plotControls.ts";

interface HeatmapConfig {
  cellSize: number;
  showValues: boolean;
}

const heatmapConfigs = new WeakMap<HTMLElement, HeatmapConfig>();

export interface HeatmapOptions {
  /** Cell values; null = empty. */
  values: (number | null)[][];
  /** Optional separate matrix of strings to display on each cell. */
  labels?: (string | null)[][];
  rowLabels: string[];
  colLabels: string[];
  /** Colour stops in [0, 1]. Default: white → orange → red ramp. */
  colorStops?: { t: number; rgb: [number, number, number] }[];
  /** Domain over which `values` are mapped to the colour ramp. */
  vmin?: number;
  vmax?: number;
  /** Caption shown above the colour bar (e.g. "−log10 p"). */
  legendTitle?: string;
  /** Cell side in px. */
  cellSize?: number;
  /** Show the numeric values inside cells (default: true). */
  showValues?: boolean;
  /** Extra CSS for cell font. */
  cellFontSize?: number;
}

const DEFAULT_STOPS = [
  { t: 0, rgb: [254, 230, 206] as [number, number, number] }, // pale peach
  { t: 0.5, rgb: [253, 174, 107] as [number, number, number] }, // orange
  { t: 1, rgb: [230, 85, 13] as [number, number, number] }, // deep orange
];

export function renderHeatmap(container: HTMLElement, opts: HeatmapOptions): void {
  // Per-container config (cell size, show-values). Initialise from the first
  // call's options so the popover reflects what the user sees.
  const cfg: HeatmapConfig = heatmapConfigs.get(container) ?? {
    cellSize: opts.cellSize ?? 32,
    showValues: opts.showValues ?? true,
  };
  heatmapConfigs.set(container, cfg);

  const fields: ControlField[] = [
    { kind: "range", key: "cellSize", label: "Cell size", min: 14, max: 80, step: 2, suffix: "px" },
    { kind: "checkbox", key: "showValues", label: "Show values" },
  ];

  const draw = (): void => renderHeatmapBody(getFigureSlot(container), opts, cfg);
  draw();
  attachToolbar(container, {
    filename: "heatmap",
    title: "Heatmap",
    fields,
    config: cfg as unknown as Record<string, unknown>,
    onChange: draw,
  });
}

function renderHeatmapBody(
  slot: HTMLElement,
  opts: HeatmapOptions,
  cfg: HeatmapConfig,
): void {
  const cell = cfg.cellSize;
  const fontSize = opts.cellFontSize ?? 10;
  const stops = opts.colorStops ?? DEFAULT_STOPS;

  // Determine value domain (ignoring null / NaN cells).
  let vmin = opts.vmin;
  let vmax = opts.vmax;
  if (vmin === undefined || vmax === undefined) {
    let lo = Infinity;
    let hi = -Infinity;
    for (const row of opts.values) {
      for (const v of row) {
        if (v === null || !Number.isFinite(v)) continue;
        if (v < lo) lo = v;
        if (v > hi) hi = v;
      }
    }
    if (!Number.isFinite(lo) || !Number.isFinite(hi)) {
      lo = 0;
      hi = 1;
    }
    if (vmin === undefined) vmin = lo;
    if (vmax === undefined) vmax = hi === lo ? lo + 1 : hi;
  }

  // Layout.
  const padTop = Math.max(60, longestStringPx(opts.colLabels, fontSize) + 12);
  const padLeft = Math.max(80, longestStringPx(opts.rowLabels, fontSize) + 12);
  const padRight = 90; // colour bar
  const padBottom = 20;
  const w = padLeft + opts.colLabels.length * cell + padRight;
  const h = padTop + opts.rowLabels.length * cell + padBottom;

  const svg: string[] = [];
  svg.push(
    `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${w} ${h}" width="${w}" height="${h}" style="font-family: Arial, sans-serif; font-size:${fontSize}px;">`,
  );

  // Column labels (rotated 90°, anchored above each column).
  for (let c = 0; c < opts.colLabels.length; c++) {
    const cx = padLeft + c * cell + cell / 2;
    const cy = padTop - 6;
    svg.push(
      `<text x="${cx}" y="${cy}" text-anchor="start" transform="rotate(-90 ${cx} ${cy})">${escapeXml(opts.colLabels[c]!)}</text>`,
    );
  }

  // Row labels.
  for (let r = 0; r < opts.rowLabels.length; r++) {
    const ry = padTop + r * cell + cell / 2 + fontSize / 3;
    svg.push(
      `<text x="${padLeft - 6}" y="${ry}" text-anchor="end">${escapeXml(opts.rowLabels[r]!)}</text>`,
    );
  }

  // Cells.
  const showValues = cfg.showValues;
  for (let r = 0; r < opts.rowLabels.length; r++) {
    for (let c = 0; c < opts.colLabels.length; c++) {
      const x = padLeft + c * cell;
      const y = padTop + r * cell;
      const v = opts.values[r]?.[c] ?? null;
      const empty = v === null || !Number.isFinite(v);
      const fill = empty ? "#f3f3f3" : rgbToHex(colorAt((v! - vmin) / (vmax - vmin), stops));
      svg.push(
        `<rect x="${x}" y="${y}" width="${cell}" height="${cell}" fill="${fill}" stroke="#fff" stroke-width="1"/>`,
      );
      if (showValues && !empty) {
        const lbl = opts.labels?.[r]?.[c] ?? formatCell(v!);
        const luma = relativeLuminance(colorAt((v! - vmin) / (vmax - vmin), stops));
        const textColor = luma > 0.6 ? "#1d1d1d" : "#ffffff";
        svg.push(
          `<text x="${x + cell / 2}" y="${y + cell / 2 + fontSize / 3}" text-anchor="middle" fill="${textColor}" style="font-size:${Math.min(fontSize, cell / 2.5)}px">${escapeXml(lbl)}</text>`,
        );
      }
    }
  }

  // Colour bar.
  const cbX = padLeft + opts.colLabels.length * cell + 30;
  const cbY = padTop;
  const cbW = 16;
  const cbH = opts.rowLabels.length * cell;
  const gradId = `g${Math.floor(Math.random() * 1e9)}`;
  svg.push(`<defs><linearGradient id="${gradId}" x1="0" y1="1" x2="0" y2="0">`);
  for (const s of stops) {
    svg.push(`<stop offset="${(s.t * 100).toFixed(0)}%" stop-color="${rgbToHex(s.rgb)}"/>`);
  }
  svg.push(`</linearGradient></defs>`);
  svg.push(
    `<rect x="${cbX}" y="${cbY}" width="${cbW}" height="${cbH}" fill="url(#${gradId})" stroke="#999"/>`,
  );
  // Bar tick marks: top, middle, bottom.
  const ticks = [vmax, (vmin + vmax) / 2, vmin];
  for (let i = 0; i < ticks.length; i++) {
    const ty = cbY + (i * cbH) / 2 + fontSize / 3;
    svg.push(
      `<text x="${cbX + cbW + 4}" y="${ty}" text-anchor="start">${formatCell(ticks[i]!)}</text>`,
    );
  }
  if (opts.legendTitle) {
    svg.push(
      `<text x="${cbX}" y="${cbY - 6}" text-anchor="start" style="font-size:${fontSize}px;">${escapeXml(opts.legendTitle)}</text>`,
    );
  }

  svg.push(`</svg>`);
  slot.innerHTML = svg.join("");
}

function colorAt(
  t: number,
  stops: { t: number; rgb: [number, number, number] }[],
): [number, number, number] {
  if (!Number.isFinite(t)) return [240, 240, 240];
  const tt = Math.max(0, Math.min(1, t));
  for (let i = 1; i < stops.length; i++) {
    const a = stops[i - 1]!;
    const b = stops[i]!;
    if (tt <= b.t) {
      const local = (tt - a.t) / (b.t - a.t || 1);
      return [
        Math.round(a.rgb[0] + (b.rgb[0] - a.rgb[0]) * local),
        Math.round(a.rgb[1] + (b.rgb[1] - a.rgb[1]) * local),
        Math.round(a.rgb[2] + (b.rgb[2] - a.rgb[2]) * local),
      ];
    }
  }
  return stops[stops.length - 1]!.rgb;
}

function rgbToHex(rgb: [number, number, number]): string {
  const h = (n: number) => n.toString(16).padStart(2, "0");
  return `#${h(rgb[0])}${h(rgb[1])}${h(rgb[2])}`;
}

function relativeLuminance(rgb: [number, number, number]): number {
  // Quick perceptual luminance — good enough for choosing dark vs light text.
  const [r, g, b] = rgb.map((c) => c / 255) as [number, number, number];
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

function formatCell(v: number): string {
  if (!Number.isFinite(v)) return "";
  if (Number.isInteger(v)) return String(v);
  if (Math.abs(v) < 0.01) return v.toExponential(1);
  if (Math.abs(v) < 1) return v.toFixed(3);
  if (Math.abs(v) < 100) return v.toFixed(2);
  return v.toFixed(0);
}

function escapeXml(s: string): string {
  return s
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&apos;");
}

function longestStringPx(strings: string[], fontSize: number): number {
  let max = 0;
  for (const s of strings) max = Math.max(max, s.length);
  // Crude but sufficient: ~0.55× font size per char in Arial.
  return Math.ceil(max * fontSize * 0.55);
}
