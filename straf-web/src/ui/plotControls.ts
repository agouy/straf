/**
 * Plot toolbar — adds a small gear + download button group to any plot
 * container. The gear opens a popover with per-plot graphical parameters;
 * the download button saves the rendered SVG (or rasterizes to PNG).
 *
 * Used by every chart helper in `plots.ts` and by the SVG heatmap.
 *
 * The toolbar lives as direct children of the container element; the actual
 * plot SVG/figure is mounted into a `.plot-figure` child so re-renders don't
 * blow away the toolbar.
 */

export type ControlField =
  | { kind: "color"; key: string; label: string }
  | { kind: "range"; key: string; label: string; min: number; max: number; step: number; suffix?: string }
  | { kind: "checkbox"; key: string; label: string }
  | {
      kind: "select";
      key: string;
      label: string;
      options: { value: string; label: string }[];
    };

export interface ToolbarOptions<T> {
  /** Stable name used in the downloaded file (e.g. "pca", "ld_heatmap"). */
  filename: string;
  /** Human-friendly title shown at the top of the popover. */
  title?: string;
  /** Configurable parameters; bound to keys on `config`. */
  fields: ControlField[];
  /** Mutable config object — updated in-place when user changes a control. */
  config: T;
  /** Called after every change so the plot can re-render. */
  onChange: () => void;
}

/**
 * Get or create the inner `.plot-figure` element. Plot helpers should mount
 * their figure here instead of replacing the entire container, so the toolbar
 * stays put.
 */
export function getFigureSlot(container: HTMLElement): HTMLElement {
  container.classList.add("plot-with-toolbar");
  let slot = container.querySelector<HTMLElement>(":scope > .plot-figure");
  if (!slot) {
    slot = document.createElement("div");
    slot.className = "plot-figure";
    container.append(slot);
  }
  return slot;
}

/**
 * Attach (or refresh) the gear + download toolbar on a plot container.
 * Idempotent: subsequent calls only update the bound config / handlers.
 */
export function attachToolbar<T extends Record<string, unknown>>(
  container: HTMLElement,
  opts: ToolbarOptions<T>,
): void {
  let toolbar = container.querySelector<HTMLElement>(":scope > .plot-toolbar");
  if (!toolbar) {
    toolbar = document.createElement("div");
    toolbar.className = "plot-toolbar";
    toolbar.innerHTML = `
      <button type="button" class="plot-tool plot-gear" title="Graphical parameters" aria-label="Graphical parameters">⚙</button>
      <button type="button" class="plot-tool plot-dl-svg" title="Download SVG" aria-label="Download SVG">SVG</button>
      <button type="button" class="plot-tool plot-dl-png" title="Download PNG" aria-label="Download PNG">PNG</button>
    `;
    // Toolbar must be the FIRST child for stable absolute positioning.
    container.prepend(toolbar);
  }

  const gearBtn = toolbar.querySelector<HTMLButtonElement>(".plot-gear")!;
  const svgBtn = toolbar.querySelector<HTMLButtonElement>(".plot-dl-svg")!;
  const pngBtn = toolbar.querySelector<HTMLButtonElement>(".plot-dl-png")!;

  // (Re)bind handlers — they capture the latest opts via closure.
  gearBtn.onclick = (e) => {
    e.stopPropagation();
    openPopover(gearBtn, opts);
  };
  svgBtn.onclick = () => downloadFigure(container, opts.filename, "svg");
  pngBtn.onclick = () => downloadFigure(container, opts.filename, "png");

  // Hide the gear if no controllable fields.
  gearBtn.style.display = opts.fields.length === 0 ? "none" : "";
}

// --- popover ---------------------------------------------------------------

let openPop: HTMLElement | null = null;

function closeOpenPopover(): void {
  if (openPop) {
    openPop.remove();
    openPop = null;
  }
}

document.addEventListener("click", (e) => {
  if (!openPop) return;
  if (openPop.contains(e.target as Node)) return;
  closeOpenPopover();
});
document.addEventListener("keydown", (e) => {
  if (e.key === "Escape") closeOpenPopover();
});

function openPopover<T extends Record<string, unknown>>(
  anchor: HTMLElement,
  opts: ToolbarOptions<T>,
): void {
  if (openPop) {
    closeOpenPopover();
    return;
  }

  const pop = document.createElement("div");
  pop.className = "plot-popover";

  const title = document.createElement("div");
  title.className = "plot-popover-title";
  title.textContent = opts.title ?? "Plot settings";
  pop.append(title);

  for (const f of opts.fields) {
    pop.append(buildField(f, opts.config, opts.onChange));
  }

  document.body.append(pop);

  // Position to the right of the anchor, bounded within the viewport.
  const rect = anchor.getBoundingClientRect();
  const popRect = pop.getBoundingClientRect();
  let left = rect.right + window.scrollX + 4;
  let top = rect.top + window.scrollY;
  if (left + popRect.width > window.innerWidth + window.scrollX - 8) {
    left = rect.left + window.scrollX - popRect.width - 4;
  }
  if (top + popRect.height > window.innerHeight + window.scrollY - 8) {
    top = window.innerHeight + window.scrollY - popRect.height - 8;
  }
  pop.style.left = `${left}px`;
  pop.style.top = `${top}px`;

  openPop = pop;
}

function buildField<T extends Record<string, unknown>>(
  field: ControlField,
  config: T,
  onChange: () => void,
): HTMLElement {
  const wrap = document.createElement("label");
  wrap.className = "plot-popover-field";

  const lbl = document.createElement("span");
  lbl.className = "plot-popover-label";
  lbl.textContent = field.label;
  wrap.append(lbl);

  switch (field.kind) {
    case "color": {
      const input = document.createElement("input");
      input.type = "color";
      input.value = String(config[field.key] ?? "#000000");
      input.addEventListener("input", () => {
        (config as Record<string, unknown>)[field.key] = input.value;
        onChange();
      });
      wrap.append(input);
      break;
    }
    case "range": {
      const row = document.createElement("div");
      row.className = "plot-popover-range";
      const input = document.createElement("input");
      input.type = "range";
      input.min = String(field.min);
      input.max = String(field.max);
      input.step = String(field.step);
      const initial = Number(config[field.key] ?? field.min);
      input.value = String(initial);
      const out = document.createElement("output");
      out.textContent = `${input.value}${field.suffix ?? ""}`;
      input.addEventListener("input", () => {
        (config as Record<string, unknown>)[field.key] = Number(input.value);
        out.textContent = `${input.value}${field.suffix ?? ""}`;
        onChange();
      });
      row.append(input, out);
      wrap.append(row);
      break;
    }
    case "checkbox": {
      const input = document.createElement("input");
      input.type = "checkbox";
      input.checked = Boolean(config[field.key]);
      input.addEventListener("change", () => {
        (config as Record<string, unknown>)[field.key] = input.checked;
        onChange();
      });
      // For checkboxes, put the input before the label text for native styling.
      wrap.classList.add("plot-popover-field-check");
      wrap.prepend(input);
      break;
    }
    case "select": {
      const sel = document.createElement("select");
      for (const o of field.options) {
        const opt = document.createElement("option");
        opt.value = o.value;
        opt.textContent = o.label;
        sel.append(opt);
      }
      sel.value = String(config[field.key] ?? field.options[0]?.value ?? "");
      sel.addEventListener("change", () => {
        (config as Record<string, unknown>)[field.key] = sel.value;
        onChange();
      });
      wrap.append(sel);
      break;
    }
  }
  return wrap;
}

// --- downloads -------------------------------------------------------------

function downloadFigure(
  container: HTMLElement,
  filename: string,
  format: "svg" | "png",
): void {
  const slot = container.querySelector<HTMLElement>(":scope > .plot-figure");
  if (!slot) {
    console.warn("Nothing to download — no .plot-figure inside container.");
    return;
  }

  const svgs = Array.from(slot.querySelectorAll<SVGSVGElement>("svg"));
  if (svgs.length === 0) {
    console.warn("Nothing to download — no SVG inside container.");
    return;
  }

  // Two cases:
  //   1. Single-SVG plot (scatter, dendrogram, heatmap, …): serialize as-is.
  //   2. Facet grid (several SVGs in a CSS grid): composite them into one SVG
  //      preserving each child's on-screen position relative to the slot.
  let svgString: string;
  let bbox: { width: number; height: number };
  if (svgs.length === 1) {
    const fig = svgs[0]!;
    svgString = serializeSvgWithStyle(fig);
    const r = fig.getBoundingClientRect();
    bbox = { width: r.width, height: r.height };
  } else {
    const composite = compositeSvgs(slot, svgs);
    svgString = composite.svgString;
    bbox = composite.bbox;
  }

  if (format === "svg") {
    triggerBlobDownload(`${filename}.svg`, svgString, "image/svg+xml;charset=utf-8");
  } else {
    rasterize(svgString, bbox).then((png) => {
      triggerBlobDownload(`${filename}.png`, png, "image/png");
    });
  }
}

/**
 * Build a single SVG that contains every child SVG of `slot` placed at its
 * on-screen position relative to the slot. Used for download of multi-figure
 * compositions (e.g. the allele frequency facet grid).
 */
function compositeSvgs(
  slot: HTMLElement,
  svgs: SVGSVGElement[],
): { svgString: string; bbox: { width: number; height: number } } {
  const slotRect = slot.getBoundingClientRect();
  const padding = 8;
  let maxX = 0;
  let maxY = 0;

  const parts: string[] = [];
  for (const fig of svgs) {
    const r = fig.getBoundingClientRect();
    const x = Math.round(r.left - slotRect.left + padding);
    const y = Math.round(r.top - slotRect.top + padding);
    const w = Math.round(r.width);
    const h = Math.round(r.height);
    maxX = Math.max(maxX, x + w);
    maxY = Math.max(maxY, y + h);

    const clone = fig.cloneNode(true) as SVGSVGElement;
    if (!clone.getAttribute("xmlns")) clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    clone.setAttribute("x", String(x));
    clone.setAttribute("y", String(y));
    clone.setAttribute("width", String(w));
    clone.setAttribute("height", String(h));
    parts.push(new XMLSerializer().serializeToString(clone));

    // Title above each sub-plot (mirror the .facet-title shown in the DOM).
    const title = fig.parentElement?.querySelector<HTMLElement>(".facet-title");
    if (title) {
      const tx = x + w / 2;
      const ty = y - 2;
      const txt = title.textContent ?? "";
      parts.unshift(
        `<text x="${tx}" y="${ty}" text-anchor="middle" font-family="Arial, sans-serif" font-size="11" font-weight="600" fill="#1d1d1d">${escapeXmlForSvg(txt)}</text>`,
      );
    }
  }

  const W = maxX + padding;
  const H = maxY + padding;
  const wrapped = [
    `<?xml version="1.0" encoding="UTF-8"?>`,
    `<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}">`,
    `<rect x="0" y="0" width="${W}" height="${H}" fill="#ffffff"/>`,
    ...parts,
    `</svg>`,
  ].join("\n");

  return { svgString: wrapped, bbox: { width: W, height: H } };
}

function escapeXmlForSvg(s: string): string {
  return s
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;");
}

/**
 * Serialize an SVG with a white background and an inline copy of computed
 * fonts/colors, so the file renders identically when reopened or rasterized.
 */
function serializeSvgWithStyle(svg: SVGSVGElement): string {
  const clone = svg.cloneNode(true) as SVGSVGElement;
  // Make sure the namespace + dimensions are explicit (Plot leaves these implicit).
  if (!clone.getAttribute("xmlns")) clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  // Inline the bbox-derived width/height so PNG export gets a sensible size.
  const bbox = svg.getBoundingClientRect();
  if (!clone.hasAttribute("width")) clone.setAttribute("width", String(Math.round(bbox.width)));
  if (!clone.hasAttribute("height")) clone.setAttribute("height", String(Math.round(bbox.height)));
  // Inject a white background rect at the bottom of the layer stack.
  const bg = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  bg.setAttribute("x", "0");
  bg.setAttribute("y", "0");
  bg.setAttribute("width", "100%");
  bg.setAttribute("height", "100%");
  bg.setAttribute("fill", "#ffffff");
  clone.insertBefore(bg, clone.firstChild);

  const xml = new XMLSerializer().serializeToString(clone);
  // Add the XML declaration so this file is portable to other tools.
  return `<?xml version="1.0" encoding="UTF-8"?>\n${xml}`;
}

function triggerBlobDownload(filename: string, data: string | Blob, mime: string): void {
  const blob = data instanceof Blob ? data : new Blob([data], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.append(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

async function rasterize(
  svgString: string,
  bbox: { width: number; height: number },
): Promise<Blob> {
  const dpr = Math.max(2, window.devicePixelRatio || 1);
  const w = Math.round(bbox.width * dpr);
  const h = Math.round(bbox.height * dpr);
  const blobUrl = URL.createObjectURL(new Blob([svgString], { type: "image/svg+xml" }));
  try {
    const img = await loadImage(blobUrl);
    const canvas = document.createElement("canvas");
    canvas.width = w;
    canvas.height = h;
    const ctx = canvas.getContext("2d")!;
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, w, h);
    ctx.drawImage(img, 0, 0, w, h);
    return await new Promise<Blob>((resolve, reject) => {
      canvas.toBlob((b) => (b ? resolve(b) : reject(new Error("toBlob failed"))), "image/png");
    });
  } finally {
    URL.revokeObjectURL(blobUrl);
  }
}

function loadImage(url: string): Promise<HTMLImageElement> {
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.onload = () => resolve(img);
    img.onerror = () => reject(new Error("Failed to load SVG image"));
    img.src = url;
  });
}
