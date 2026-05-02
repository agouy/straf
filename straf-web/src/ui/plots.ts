/**
 * Plot helpers built on Observable Plot.
 *
 * Each helper:
 *   1. mounts the figure into a `.plot-figure` slot inside the container
 *      (so the gear-button toolbar isn't wiped on re-render),
 *   2. maintains a per-container config object whose values can be edited
 *      via the toolbar's gear popover, and
 *   3. wires SVG/PNG download buttons.
 *
 * Public API matches the original Plotly-backed helpers so call sites in
 * main.ts don't need to change.
 */

import * as Plot from "@observablehq/plot";
import { attachToolbar, getFigureSlot } from "./plotControls.ts";
import type { ControlField } from "./plotControls.ts";

const PALETTE = [
  "#dd4814", "#3b8ea5", "#5cab7d", "#a86ec9", "#e2a83a",
  "#666666", "#c95792", "#3a6cd2", "#94a800", "#b04848",
  "#3aa2c0", "#7d6f53", "#bf6b50", "#5fa8d3", "#ce9b5f",
];

export function paletteColor(i: number): string {
  return PALETTE[i % PALETTE.length]!;
}

// Each plot type maintains its own config map keyed by container element.
const barConfigs = new WeakMap<HTMLElement, BarConfig>();
const gridConfigs = new WeakMap<HTMLElement, GridConfig>();
const scatterConfigs = new WeakMap<HTMLElement, ScatterConfig>();
const labeledConfigs = new WeakMap<HTMLElement, ScatterConfig>();
const dendroConfigs = new WeakMap<HTMLElement, DendroConfig>();

interface BarConfig {
  color: string;
  width: number;
  height: number;
  showValues: boolean;
}
interface GridConfig {
  color: string;
  width: number;
  height: number;
}
interface ScatterConfig {
  pointSize: number;
  width: number;
  height: number;
  showLabels: boolean;
}
interface DendroConfig {
  width: number;
  lineWidth: number;
  fontSize: number;
}

// --- bar plot --------------------------------------------------------------

export function barPlot(
  container: HTMLElement,
  labels: string[],
  values: number[],
  xlabel: string,
  color = "#36648B",
): void {
  const cfg: BarConfig = barConfigs.get(container) ?? {
    color,
    width: 720,
    height: Math.max(300, 22 * labels.length + 80),
    showValues: false,
  };
  barConfigs.set(container, cfg);

  const fields: ControlField[] = [
    { kind: "color", key: "color", label: "Bar color" },
    { kind: "range", key: "width", label: "Plot width", min: 400, max: 1400, step: 20, suffix: "px" },
    { kind: "range", key: "height", label: "Plot height", min: 200, max: 1200, step: 20, suffix: "px" },
    { kind: "checkbox", key: "showValues", label: "Show values on bars" },
  ];

  const render = (): void => {
    const data = labels.map((label, i) => ({ label, value: values[i] ?? 0 }));
    const fig = Plot.plot({
      marginLeft: Math.max(110, longestStringPx(labels) + 24),
      marginBottom: 40,
      marginRight: 20,
      marginTop: 20,
      width: cfg.width,
      height: cfg.height,
      x: { label: xlabel, grid: true },
      y: { label: null, domain: labels },
      marks: [
        Plot.barX(data, {
          y: "label",
          x: "value",
          fill: cfg.color,
          tip: true,
          title: (d: { label: string; value: number }) =>
            `${d.label}\n${xlabel}: ${formatNum(d.value)}`,
        }),
        ...(cfg.showValues
          ? [
              Plot.text(data, {
                y: "label",
                x: "value",
                text: (d: { value: number }) => formatNum(d.value),
                dx: 4,
                textAnchor: "start",
                fontSize: 10,
                fill: "#333",
              }),
            ]
          : []),
        Plot.ruleX([0]),
      ],
    });
    getFigureSlot(container).replaceChildren(fig);
  };

  render();
  attachToolbar(container, {
    filename: "barplot",
    title: "Bar plot",
    fields,
    config: cfg as unknown as Record<string, unknown>,
    onChange: render,
  });
}

// --- multi-locus allele frequency grid ------------------------------------

export function alleleFreqGrid(
  container: HTMLElement,
  loci: string[],
  freqs: { allele: string; count: number }[][],
  color = "#36648B",
): void {
  const cfg: GridConfig = gridConfigs.get(container) ?? {
    color,
    width: 1100,
    height: Math.max(240, Math.ceil(loci.length / (2 + Math.floor(loci.length / 5))) * 180),
  };
  gridConfigs.set(container, cfg);

  const fields: ControlField[] = [
    { kind: "color", key: "color", label: "Bar color" },
    { kind: "range", key: "width", label: "Plot width", min: 600, max: 1600, step: 20, suffix: "px" },
    { kind: "range", key: "height", label: "Plot height", min: 200, max: 1600, step: 40, suffix: "px" },
  ];

  const render = (): void => {
    // Plot's `fx` faceting shares the x-axis across facets, which crushes the
    // bars when each locus has its own allele set. We render one small Plot
    // per locus instead and arrange them in a CSS grid — same effect as
    // ggplot's `facet_wrap(scales = "free_x")`.
    const nL = loci.length;
    const nCol = Math.min(nL, Math.max(2, Math.round(cfg.width / 240)));
    const nRow = Math.ceil(nL / nCol);
    const subWidth = Math.floor(cfg.width / nCol);
    const subHeight = Math.max(140, Math.floor(cfg.height / nRow));

    const grid = document.createElement("div");
    grid.style.display = "grid";
    grid.style.gridTemplateColumns = `repeat(${nCol}, 1fr)`;
    grid.style.gap = "4px";
    grid.style.width = `${cfg.width}px`;

    for (let i = 0; i < nL; i++) {
      const locus = loci[i]!;
      const cell = document.createElement("div");
      cell.className = "facet-cell";

      const title = document.createElement("div");
      title.className = "facet-title";
      title.textContent = locus;
      cell.append(title);

      const sub = Plot.plot({
        width: subWidth,
        height: subHeight,
        marginLeft: 36,
        marginRight: 6,
        marginTop: 6,
        marginBottom: 28,
        x: { label: null, type: "band", padding: 0.15, tickRotate: -45 },
        y: { label: null, grid: true },
        marks: [
          Plot.barY(freqs[i]!, {
            x: "allele",
            y: "count",
            fill: cfg.color,
            tip: true,
            title: (d: { allele: string; count: number }) =>
              `${locus} · ${d.allele}\nn = ${formatNum(d.count)}`,
          }),
          Plot.ruleY([0]),
        ],
      });
      cell.append(sub);
      grid.append(cell);
    }
    getFigureSlot(container).replaceChildren(grid);
  };

  render();
  attachToolbar(container, {
    filename: "allele_frequency_grid",
    title: "Allele frequency grid",
    fields,
    config: cfg as unknown as Record<string, unknown>,
    onChange: render,
  });
}

// --- scatter (PCA-style: groups + ellipses) -------------------------------

export interface ScatterGroup {
  name: string;
  x: number[];
  y: number[];
  text: string[];
}

export function scatterPlot(
  container: HTMLElement,
  groups: ScatterGroup[],
  xlabel: string,
  ylabel: string,
  ellipses?: { name: string; x: number[]; y: number[] }[],
): void {
  const cfg: ScatterConfig = scatterConfigs.get(container) ?? {
    pointSize: 4,
    width: 720,
    height: 480,
    showLabels: false,
  };
  scatterConfigs.set(container, cfg);

  const fields: ControlField[] = [
    { kind: "range", key: "pointSize", label: "Point size", min: 2, max: 12, step: 1, suffix: "px" },
    { kind: "range", key: "width", label: "Plot width", min: 400, max: 1200, step: 20, suffix: "px" },
    { kind: "range", key: "height", label: "Plot height", min: 300, max: 900, step: 20, suffix: "px" },
    { kind: "checkbox", key: "showLabels", label: "Show point labels" },
  ];

  const render = (): void => {
    const points: { group: string; x: number; y: number; label: string }[] = [];
    for (const g of groups) {
      for (let i = 0; i < g.x.length; i++) {
        points.push({ group: g.name, x: g.x[i]!, y: g.y[i]!, label: g.text[i] ?? "" });
      }
    }
    const ellipseData: { group: string; x: number; y: number }[] = [];
    if (ellipses) {
      for (const e of ellipses) {
        for (let i = 0; i < e.x.length; i++) {
          ellipseData.push({ group: e.name, x: e.x[i]!, y: e.y[i]! });
        }
      }
    }
    const groupNames = groups.map((g) => g.name);
    const colorRange = groupNames.map((_, i) => paletteColor(i));

    const fig = Plot.plot({
      width: cfg.width,
      height: cfg.height,
      marginLeft: 60,
      marginBottom: 50,
      marginRight: 20,
      marginTop: 20,
      x: { label: xlabel, grid: true },
      y: { label: ylabel, grid: true },
      color: { domain: groupNames, range: colorRange, legend: true },
      marks: [
        Plot.ruleX([0], { stroke: "#ccc" }),
        Plot.ruleY([0], { stroke: "#ccc" }),
        ...(ellipseData.length > 0
          ? [
              Plot.line(ellipseData, {
                x: "x",
                y: "y",
                stroke: "group",
                z: "group",
                strokeWidth: 1.2,
                strokeOpacity: 0.9,
              }),
            ]
          : []),
        Plot.dot(points, {
          x: "x",
          y: "y",
          fill: "group",
          stroke: "#222",
          strokeWidth: 0.4,
          r: cfg.pointSize,
          tip: true,
          title: (d: { label: string; group: string; x: number; y: number }) =>
            `${d.label}\n${d.group}\n${xlabel}: ${formatNum(d.x)}\n${ylabel}: ${formatNum(d.y)}`,
        }),
        ...(cfg.showLabels
          ? [
              Plot.text(points, {
                x: "x",
                y: "y",
                text: "label",
                dy: -10,
                fontSize: 10,
                fill: "group",
              }),
            ]
          : []),
      ],
    });
    getFigureSlot(container).replaceChildren(fig);
  };

  render();
  attachToolbar(container, {
    filename: "scatter",
    title: "Scatter plot",
    fields,
    config: cfg as unknown as Record<string, unknown>,
    onChange: render,
  });
}

// --- labeled scatter (one labeled point per row, e.g. populations) --------

export function labeledScatter(
  container: HTMLElement,
  labels: string[],
  xs: number[],
  ys: number[],
  xlabel: string,
  ylabel: string,
): void {
  const cfg: ScatterConfig = labeledConfigs.get(container) ?? {
    pointSize: 6,
    width: 720,
    height: 480,
    showLabels: true,
  };
  labeledConfigs.set(container, cfg);

  const fields: ControlField[] = [
    { kind: "range", key: "pointSize", label: "Point size", min: 3, max: 14, step: 1, suffix: "px" },
    { kind: "range", key: "width", label: "Plot width", min: 400, max: 1200, step: 20, suffix: "px" },
    { kind: "range", key: "height", label: "Plot height", min: 300, max: 900, step: 20, suffix: "px" },
    { kind: "checkbox", key: "showLabels", label: "Show labels" },
  ];

  const render = (): void => {
    const data = labels.map((label, i) => ({
      label,
      x: xs[i] ?? 0,
      y: ys[i] ?? 0,
      color: paletteColor(i),
    }));
    const fig = Plot.plot({
      width: cfg.width,
      height: cfg.height,
      marginLeft: 60,
      marginBottom: 50,
      marginRight: 20,
      marginTop: 20,
      x: { label: xlabel, grid: true },
      y: { label: ylabel, grid: true },
      color: { domain: labels, range: data.map((d) => d.color) },
      marks: [
        Plot.ruleX([0], { stroke: "#ccc" }),
        Plot.ruleY([0], { stroke: "#ccc" }),
        Plot.dot(data, {
          x: "x",
          y: "y",
          fill: "label",
          stroke: "#222",
          strokeWidth: 0.5,
          r: cfg.pointSize,
          tip: true,
          title: (d: { label: string; x: number; y: number }) =>
            `${d.label}\n${xlabel}: ${formatNum(d.x)}\n${ylabel}: ${formatNum(d.y)}`,
        }),
        ...(cfg.showLabels
          ? [
              Plot.text(data, {
                x: "x",
                y: "y",
                text: "label",
                dy: -10,
                fontSize: 11,
              }),
            ]
          : []),
      ],
    });
    getFigureSlot(container).replaceChildren(fig);
  };

  render();
  attachToolbar(container, {
    filename: "scatter",
    title: "Scatter plot",
    fields,
    config: cfg as unknown as Record<string, unknown>,
    onChange: render,
  });
}

// --- horizontal dendrogram (root left, leaves right) ----------------------

/**
 * Render a dendrogram horizontally — root on the left, leaves on the right.
 * Input segments come from `dendrogramSegments` which produces U-shaped joins
 * in (x = leaf-position, y = distance) coordinates with leaves at y = 0.
 *
 * To rotate to horizontal we swap axes: leaves end up at x = `height`, root
 * progresses leftward as merges accumulate. Segments are mirrored so the tree
 * grows naturally to the right.
 */
export function dendrogramPlot(
  container: HTMLElement,
  leafX: number[],
  leafLabels: string[],
  segments: { x: number[]; y: number[] }[],
  totalHeight: number,
): void {
  const cfg: DendroConfig = dendroConfigs.get(container) ?? {
    width: 720,
    lineWidth: 1.5,
    fontSize: 12,
  };
  dendroConfigs.set(container, cfg);

  const fields: ControlField[] = [
    { kind: "range", key: "width", label: "Plot width", min: 400, max: 1200, step: 20, suffix: "px" },
    { kind: "range", key: "lineWidth", label: "Line thickness", min: 1, max: 4, step: 0.5, suffix: "px" },
    { kind: "range", key: "fontSize", label: "Label size", min: 8, max: 18, step: 1, suffix: "px" },
  ];

  const render = (): void => {
    // Swap x/y so the tree is horizontal. Root is at x = totalHeight,
    // leaves at x = 0; we'll flip the x scale to draw root → leaves left → right.
    const lineData: { seg: number; x: number; y: number }[] = [];
    for (let s = 0; s < segments.length; s++) {
      const segment = segments[s]!;
      for (let i = 0; i < segment.x.length; i++) {
        lineData.push({
          seg: s,
          x: segment.y[i]!, // distance → horizontal
          y: segment.x[i]!, // leaf position → vertical
        });
      }
    }
    const leafData = leafLabels.map((label, i) => ({
      label,
      y: leafX[i] ?? i,
    }));

    const height = Math.max(220, leafLabels.length * (cfg.fontSize + 6) + 60);
    const fig = Plot.plot({
      width: cfg.width,
      height,
      marginLeft: 30,
      marginRight: Math.max(80, longestStringPx(leafLabels, cfg.fontSize) + 20),
      marginTop: 20,
      marginBottom: 40,
      x: {
        label: "Distance",
        domain: [totalHeight * 1.02, 0], // reversed: root left, leaves right
        grid: true,
      },
      y: {
        label: null,
        domain: [-0.5, leafLabels.length - 0.5],
        ticks: [],
      },
      marks: [
        Plot.line(lineData, {
          x: "x",
          y: "y",
          z: "seg",
          stroke: "#444",
          strokeWidth: cfg.lineWidth,
        }),
        Plot.dot(leafData, {
          x: 0,
          y: "y",
          r: 2.5,
          fill: "#dd4814",
          tip: true,
          title: (d: { label: string }) => d.label,
        }),
        Plot.text(leafData, {
          x: 0,
          y: "y",
          text: "label",
          dx: 6,
          textAnchor: "start",
          fontSize: cfg.fontSize,
        }),
      ],
    });
    getFigureSlot(container).replaceChildren(fig);
  };

  render();
  attachToolbar(container, {
    filename: "dendrogram",
    title: "Dendrogram",
    fields,
    config: cfg as unknown as Record<string, unknown>,
    onChange: render,
  });
}

// --- helpers ---------------------------------------------------------------

function formatNum(v: number): string {
  if (!Number.isFinite(v)) return "";
  if (Number.isInteger(v)) return String(v);
  if (Math.abs(v) < 0.01) return v.toExponential(2);
  if (Math.abs(v) < 1) return v.toFixed(4);
  if (Math.abs(v) < 100) return v.toFixed(3);
  return v.toFixed(1);
}

function longestStringPx(strings: string[], fontSize = 11): number {
  let max = 0;
  for (const s of strings) max = Math.max(max, s.length);
  return Math.ceil(max * fontSize * 0.55);
}
