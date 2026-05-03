import {
  parseStrafTsv,
  allMask,
  popMask,
  uniquePopulations,
  ParseError,
  getParseWarnings,
} from "./parser.ts";
import type { Genotypes, Ploidy, ParseIssue } from "./parser.ts";
import { alleleFrequencies, perLocusFrequencies } from "./stats/freq.ts";
import { locusIndices, combinedIndices } from "./stats/forensic.ts";
import type { LocusStats } from "./stats/forensic.ts";
import { hweChiSquare } from "./stats/hwe.ts";
import { pca } from "./stats/pca.ts";
import type { PcaResult } from "./stats/pca.ts";
import { populationDistances } from "./stats/distances.ts";
import type { DistanceMethod } from "./stats/distances.ts";
import { classicalMds } from "./stats/mds.ts";
import { confidenceEllipse } from "./stats/ellipse.ts";
import { pairwiseFst, perLocusFStats, basicStatsPerLocus } from "./stats/fst.ts";
import type { LocusFStats, LocusBasicStats } from "./stats/fst.ts";
import { pairwiseLd, ldMatrix } from "./stats/ld.ts";
import { isoMds } from "./stats/iso_mds.ts";
import { haplotypeStatsForPop } from "./stats/haplotype.ts";
import { writeXlsx } from "./ui/xlsx.ts";
import { renderHeatmap } from "./ui/heatmap.ts";
import {
  parseReferenceCsv,
  dropMissingColumns,
  neiDistanceFromMatrix,
} from "./stats/reference.ts";
import type { ReferenceFrequencies } from "./stats/reference.ts";
import { renderTable, fmtNum, downloadText, tsvFromMatrix } from "./ui/table.ts";
import {
  barPlot,
  alleleFreqGrid,
  scatterPlot,
  labeledScatter,
  dendrogramPlot,
} from "./ui/plots.ts";
import { hclust, dendrogramSegments } from "./stats/hclust.ts";
import {
  toGenepop,
  toArlequin,
  toFamilias,
  toAlleleFreqCsv,
} from "./convert/formats.ts";

interface AppState {
  genos: Genotypes | null;
  pcaResult: PcaResult | null;
  rawText: string | null;
  /** Lazy sections the user has triggered — replayed when the dataset is reloaded. */
  ran: { pca: boolean; mds: boolean; ld: boolean; fst: boolean };
}

const state: AppState = {
  genos: null,
  pcaResult: null,
  rawText: null,
  ran: { pca: false, mds: false, ld: false, fst: false },
};

/**
 * Defer heavy compute until after the next browser paint, so any UI placed
 * before the call (spinner, status text) actually becomes visible to the user.
 * A single rAF runs *before* the upcoming paint; the second rAF runs after it.
 */
function deferCompute(fn: () => void): void {
  requestAnimationFrame(() => requestAnimationFrame(fn));
}

const $ = <T extends HTMLElement = HTMLElement>(id: string): T =>
  document.getElementById(id) as T;

// --- Top-level page tabs ---------------------------------------------------
document.querySelectorAll<HTMLButtonElement>(".navbar-tab").forEach((btn) => {
  btn.addEventListener("click", () => {
    const target = btn.dataset.toptab!;
    document.querySelectorAll(".navbar-tab").forEach((b) => b.classList.remove("active"));
    document.querySelectorAll(".page").forEach((p) => p.classList.remove("active"));
    btn.classList.add("active");
    $(`page-${target}`).classList.add("active");
  });
});

// --- Sub-tabs --------------------------------------------------------------
document.querySelectorAll<HTMLButtonElement>(".subtab").forEach((btn) => {
  btn.addEventListener("click", () => {
    const target = btn.dataset.subtab!;
    document.querySelectorAll(".subtab").forEach((b) => b.classList.remove("active"));
    document.querySelectorAll(".subpanel").forEach((p) => p.classList.remove("active"));
    btn.classList.add("active");
    $(`sub-${target}`).classList.add("active");
  });
});

// --- File loading ----------------------------------------------------------
$<HTMLInputElement>("file").addEventListener("change", async (e) => {
  const input = e.target as HTMLInputElement;
  const file = input.files?.[0];
  if (!file) return;
  const text = await file.text();
  loadText(text, file.name);
});

$<HTMLButtonElement>("loadExample").addEventListener("click", async () => {
  setStatus("Loading example…");
  try {
    // Pick example based on currently-selected ploidy
    const ploidy = currentPloidy();
    const fname = ploidy === 1 ? "exampleSTRAFhaplo.txt" : "exampleSTRAFdiplo.txt";
    const res = await fetch(`./${fname}`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const text = await res.text();
    loadText(text, fname);
  } catch (err) {
    setStatus(`Could not load example: ${(err as Error).message}`, "error");
  }
});

// Re-parse when ploidy changes (if data already loaded).
document.querySelectorAll<HTMLInputElement>('input[name="ploidy"]').forEach((r) => {
  r.addEventListener("change", () => {
    if (state.rawText) loadText(state.rawText, "(re-parsed)");
  });
});

function currentPloidy(): Ploidy {
  const r = document.querySelector<HTMLInputElement>('input[name="ploidy"]:checked');
  return Number(r?.value ?? 2) as Ploidy;
}

function loadText(text: string, fname: string): void {
  const ploidy = currentPloidy();
  try {
    const genos = parseStrafTsv(text, ploidy);
    state.genos = genos;
    state.rawText = text;
    state.pcaResult = null;
    setStatus(
      `Loaded ${fname}: ${genos.individuals.length} individuals, ${genos.loci.length} loci, ${uniquePopulations(genos).length} population(s).`,
      "ok",
    );
    $("noData").style.display = "none";
    $("hasData").style.display = "block";
    renderIssues(getParseWarnings(genos));
    refreshAll();
  } catch (err) {
    state.genos = null;
    $("hasData").style.display = "none";
    $("noData").style.display = "block";
    if (err instanceof ParseError) {
      setStatus(`Could not load ${fname}: ${err.issues.length} issue(s).`, "error");
      renderIssues(err.issues);
    } else {
      setStatus(`Unexpected error: ${(err as Error).message}`, "error");
      renderIssues([{ kind: "error", message: (err as Error).message }]);
    }
  }
}

function renderIssues(issues: ParseIssue[]): void {
  const el = $("parseIssues");
  if (issues.length === 0) {
    el.style.display = "none";
    el.innerHTML = "";
    return;
  }
  const hasErrors = issues.some((i) => i.kind === "error");
  el.className = "issues" + (hasErrors ? " has-errors" : "");
  el.style.display = "";
  const errCount = issues.filter((i) => i.kind === "error").length;
  const warnCount = issues.length - errCount;
  const title =
    errCount > 0
      ? `${errCount} error${errCount === 1 ? "" : "s"}` +
        (warnCount > 0 ? ` and ${warnCount} warning${warnCount === 1 ? "" : "s"}` : "")
      : `${warnCount} warning${warnCount === 1 ? "" : "s"}`;

  const groups = groupIssues(issues);
  const items = groups.map(renderIssueGroup).join("");

  el.innerHTML = `
    <button class="dismiss" type="button" aria-label="Dismiss">×</button>
    <h4>Input file: ${title}</h4>
    <ul>${items}</ul>
  `;
  const dismiss = el.querySelector<HTMLButtonElement>(".dismiss");
  if (dismiss) dismiss.addEventListener("click", () => renderIssues([]));
}

interface IssueGroup {
  kind: "error" | "warning";
  template: string;
  hint?: string;
  members: ParseIssue[];
  /** Distinct quoted values extracted from each member's message (e.g. individual IDs). */
  values: string[];
}

/** Replace quoted values in a message with "…" so messages that differ only by ID collapse together. */
function templateOf(msg: string): string {
  return msg.replace(/"[^"]*"/g, '"…"');
}

function extractQuoted(msg: string): string[] {
  const out: string[] = [];
  const re = /"([^"]*)"/g;
  let m: RegExpExecArray | null;
  while ((m = re.exec(msg)) !== null) out.push(m[1]!);
  return out;
}

function groupIssues(issues: ParseIssue[]): IssueGroup[] {
  const map = new Map<string, IssueGroup>();
  const order: string[] = [];
  for (const i of issues) {
    const tmpl = templateOf(i.message);
    const key = `${i.kind}|${tmpl}|${i.hint ?? ""}`;
    let g = map.get(key);
    if (!g) {
      g = { kind: i.kind, template: tmpl, hint: i.hint, members: [], values: [] };
      map.set(key, g);
      order.push(key);
    }
    g.members.push(i);
    for (const v of extractQuoted(i.message)) g.values.push(v);
  }
  return order.map((k) => map.get(k)!);
}

function renderIssueGroup(g: IssueGroup): string {
  const first = g.members[0]!;
  const where =
    first.line !== undefined
      ? `<span class="where">L${first.line}${first.column !== undefined ? `:C${first.column}` : ""}</span>`
      : "";
  const hint = g.hint ? `<br/><span class="hint">${escapeHtml(g.hint)}</span>` : "";

  if (g.members.length === 1) {
    return `<li class="${g.kind}">${where}${escapeHtml(first.message)}${hint}</li>`;
  }

  const lastLine = g.members[g.members.length - 1]!.line;
  const lineRange =
    first.line !== undefined && lastLine !== undefined && lastLine !== first.line
      ? `lines ${first.line}–${lastLine}`
      : first.line !== undefined
        ? `line ${first.line}`
        : "";
  const occ = `${g.members.length} occurrences${lineRange ? ` (${lineRange})` : ""}`;

  // Show up to 5 distinct example values, if the message had any quoted bits.
  const distinct: string[] = [];
  const seen = new Set<string>();
  for (const v of g.values) {
    if (seen.has(v)) continue;
    seen.add(v);
    distinct.push(v);
    if (distinct.length >= 5) break;
  }
  let examples = "";
  if (distinct.length > 0) {
    const totalDistinct = new Set(g.values).size;
    const remaining = totalDistinct - distinct.length;
    const list = distinct.map((v) => `"${escapeHtml(v)}"`).join(", ");
    examples = `<br/><span class="hint">Examples: ${list}${remaining > 0 ? ` (+${remaining} more)` : ""}</span>`;
  }

  return `<li class="${g.kind}">${where}<span class="muted small">[${occ}]</span> ${escapeHtml(g.template)}${examples}${hint}</li>`;
}

function escapeHtml(s: string): string {
  return s
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&#39;");
}

function setStatus(msg: string, kind: "ok" | "error" | "" = ""): void {
  const el = $("status");
  el.textContent = msg;
  el.className = "status" + (kind ? ` ${kind}` : "");
}

/** Show or hide a spinner badge inside a wrap element. Pass null to remove. */
function setComputing(wrapId: string, label: string | null): void {
  const wrap = $(wrapId);
  let s = wrap.querySelector<HTMLElement>(":scope > .computing");
  if (label === null) {
    if (s) s.remove();
    return;
  }
  if (!s) {
    s = document.createElement("div");
    s.className = "computing";
    wrap.prepend(s);
  }
  s.textContent = label;
}

// --- Refresh orchestrator --------------------------------------------------
function refreshAll(): void {
  if (!state.genos) return;
  populatePopSelectors();
  renderDataSummary();
  renderDataPreview();
  renderAlleleFreqPlot();
  renderFreqTable();
  renderForensic();
  renderPopgen();
  renderHaploSection();
  // PCA / MDS / LD / Fst are gated behind run buttons (lazy compute).
  // If the user had already run them on the previous dataset, replay them on
  // the new one so they don't have to re-click after every upload.
  state.pcaResult = null;
  cachedMds = null;
  syncLazySection("pca", "pcaWrap", "runPCA", computePca);
  syncLazySection("mds", "mdsWrap", "runMDS", computeMds);
  syncLazySection("ld", "ldWrap", "runLD", renderLd);
  syncLazySection("fst", "fstWrap", "runFst", renderFst);
}

function syncLazySection(
  key: keyof AppState["ran"],
  wrapId: string,
  buttonId: string,
  run: () => void,
): void {
  const wrap = $(wrapId);
  const btn = $<HTMLButtonElement>(buttonId);
  if (state.ran[key]) {
    wrap.style.display = "";
    btn.disabled = true;
    run();
  } else {
    wrap.style.display = "none";
    btn.disabled = false;
  }
}

function populatePopSelectors(): void {
  if (!state.genos) return;
  const pops = ["all", ...uniquePopulations(state.genos)];
  for (const id of ["freqPop", "forensicPop", "popgenPop", "convPop"]) {
    const sel = $<HTMLSelectElement>(id);
    if (!sel) continue;
    sel.innerHTML = pops.map((p) => `<option value="${p}">${p}</option>`).join("");
  }
}

function maskFor(pop: string): boolean[] {
  if (!state.genos) return [];
  return pop === "all" ? allMask(state.genos) : popMask(state.genos, pop);
}

// --- DATA tab --------------------------------------------------------------
$<HTMLSelectElement>("freqPop").addEventListener("change", renderFreqTable);

function renderDataSummary(): void {
  if (!state.genos) return;
  const g = state.genos;
  const pops = uniquePopulations(g);
  // Count missing genotypes (any allele null at any locus per individual).
  let missing = 0;
  let totalGeno = 0;
  for (let l = 0; l < g.loci.length; l++) {
    for (const tuple of g.alleles[l]!) {
      totalGeno++;
      if (tuple.some((a) => a === null)) missing++;
    }
  }
  const missingPct = totalGeno > 0 ? (100 * missing) / totalGeno : 0;
  // Total distinct alleles across all loci.
  let totalAlleles = 0;
  for (let l = 0; l < g.loci.length; l++) {
    const seen = new Set<string>();
    for (const tuple of g.alleles[l]!) for (const a of tuple) if (a !== null) seen.add(a);
    totalAlleles += seen.size;
  }

  const cards: { label: string; value: string; sub?: string }[] = [
    { label: "Individuals", value: String(g.individuals.length) },
    { label: "Loci", value: String(g.loci.length) },
    { label: "Populations", value: String(pops.length) },
    { label: "Ploidy", value: g.ploidy === 2 ? "Diploid" : "Haploid" },
    { label: "Distinct alleles", value: String(totalAlleles), sub: "summed across loci" },
    {
      label: "Missing genotypes",
      value: `${missing}`,
      sub: `${missingPct.toFixed(1)}% of ${totalGeno}`,
    },
  ];
  $("summaryCards").innerHTML = cards
    .map(
      (c) => `
        <div class="summary-card">
          <div class="label">${c.label}</div>
          <div class="value">${c.value}</div>
          ${c.sub ? `<div class="sub">${c.sub}</div>` : ""}
        </div>`,
    )
    .join("");

  // Per-population sizes table (rows = pops, cols = N indiv, % of total).
  const headers = ["Population", "N", "% of total"];
  const total = g.individuals.length;
  const rows: (string | number | null)[][] = pops.map((p) => {
    const n = g.populations.filter((q) => q === p).length;
    return [p, n, (100 * n) / total];
  });
  rows.push(["Total", total, 100]);
  renderTable($("popSizesTable"), headers, rows, (v, col) => {
    if (col === 0) return v as string;
    if (typeof v !== "number") return "";
    if (col === 1) return String(v);
    return v.toFixed(1) + "%";
  });
}

function renderDataPreview(): void {
  if (!state.genos) return;
  const g = state.genos;
  // Show the raw header + data as a table — matches original's DT::dataTableOutput.
  const headers = ["ind", "pop"];
  for (const l of g.loci) {
    for (let k = 0; k < g.ploidy; k++) headers.push(l);
  }
  const rows: (string | number | null)[][] = g.individuals.map((ind, i) => {
    const row: (string | number | null)[] = [ind, g.populations[i]!];
    for (let l = 0; l < g.loci.length; l++) {
      for (const a of g.alleles[l]![i]!) row.push(a ?? "0");
    }
    return row;
  });
  renderTable($("dataPreview"), headers, rows, (v) => (v === null ? "" : String(v)));
}

function renderAlleleFreqPlot(): void {
  if (!state.genos) return;
  const { freqs, N } = perLocusFrequencies(state.genos);
  const freqArr = freqs.map((m, l) => {
    const Nl = N[l]!;
    const arr = Array.from(m.entries())
      .map(([allele, p]) => ({ allele, count: Nl > 0 ? p * Nl : 0 }))
      .sort((a, b) => {
        const na = Number(a.allele);
        const nb = Number(b.allele);
        if (Number.isFinite(na) && Number.isFinite(nb)) return na - nb;
        return a.allele.localeCompare(b.allele);
      });
    return arr;
  });
  alleleFreqGrid($("alleleFreqPlot"), state.genos.loci, freqArr);
}

function renderFreqTable(): void {
  if (!state.genos) return;
  const pop = $<HTMLSelectElement>("freqPop").value || "all";
  const tbl = alleleFrequencies(state.genos, maskFor(pop));
  const headers = ["Allele", ...tbl.loci];
  const rows: (string | number | null)[][] = tbl.alleles.map((a, ai) => [
    a,
    ...tbl.values[ai]!.map((v) => v),
  ]);
  rows.push(["N", ...tbl.N]);

  renderTable($("freqTable"), headers, rows, (v, col) => {
    if (col === 0) return v as string;
    if (v === null || v === undefined) return "";
    if (typeof v === "number") {
      return Number.isInteger(v) ? String(v) : (v as number).toFixed(3);
    }
    return String(v);
  });

  $<HTMLButtonElement>("dlFreq").onclick = () => {
    downloadText(
      `allele_frequencies_${pop}.tsv`,
      tsvFromMatrix(headers, rows),
      "text/tab-separated-values",
    );
  };
  $<HTMLButtonElement>("dlFreqXL").onclick = () => {
    writeXlsx(`allele_frequencies_${pop}.xlsx`, "allele_frequencies", headers, rows);
  };
}

// --- FORENSIC tab ----------------------------------------------------------
$<HTMLInputElement>("displayForensics").addEventListener("change", (e) => {
  $("forensicWrap").style.display = (e.target as HTMLInputElement).checked ? "" : "none";
  if ((e.target as HTMLInputElement).checked) renderForensic();
});
$<HTMLSelectElement>("forensicPop").addEventListener("change", renderForensic);
$<HTMLSelectElement>("forensicCol").addEventListener("change", () => renderForensicPlot());

let lastForensicStats: LocusStats[] = [];

function renderForensic(): void {
  if (!state.genos) return;
  if (!$<HTMLInputElement>("displayForensics").checked) {
    $("forensicWrap").style.display = "none";
    return;
  }
  $("forensicWrap").style.display = "";

  const pop = $<HTMLSelectElement>("forensicPop").value || "all";
  const stats = locusIndices(state.genos, maskFor(pop));
  lastForensicStats = stats;

  // Forensic columns only — popgen-only stats live on the popgen tab.
  type Col = { key: keyof LocusStats; label: string; digits: number };
  const cols: Col[] = [
    { key: "N", label: "N", digits: 0 },
    { key: "Nall", label: "Nall", digits: 0 },
    { key: "GD", label: "GD (Hexp)", digits: 4 },
    { key: "PIC", label: "PIC", digits: 4 },
    { key: "PM", label: "PM", digits: 4 },
    { key: "PD", label: "PD", digits: 4 },
  ];
  if (state.genos.ploidy === 2) {
    cols.push(
      { key: "PE", label: "PE", digits: 4 },
      { key: "TPI", label: "TPI", digits: 4 },
    );
  }

  const headers = ["Locus", ...cols.map((c) => c.label)];
  const rows: (string | number | null)[][] = stats.map((s) => [
    s.locus,
    ...cols.map((c) => {
      const v = s[c.key];
      return typeof v === "number" || v === null || v === undefined ? (v ?? null) : null;
    }),
  ]);

  // Combined row.
  const combined = combinedIndices(stats);
  const combinedRow: (string | number | null)[] = ["Combined"];
  for (const c of cols) {
    if (c.key === "PM") combinedRow.push(combined.combinedPM);
    else if (c.key === "PD") combinedRow.push(combined.combinedPD);
    else if (c.key === "PE") combinedRow.push(combined.combinedPE ?? null);
    else combinedRow.push(null);
  }
  rows.push(combinedRow);

  renderTable($("forensicTable"), headers, rows, (v, col) => {
    if (col === 0) return v as string;
    if (v === null || v === undefined) return "";
    if (typeof v !== "number") return String(v);
    return fmtNum(v, cols[col - 1]!.digits);
  });

  // Plot column selector.
  const sel = $<HTMLSelectElement>("forensicCol");
  const prev = sel.value;
  sel.innerHTML = cols.map((c) => `<option value="${c.key}">${c.label}</option>`).join("");
  if (cols.some((c) => c.key === prev)) sel.value = prev;
  else sel.value = "Nall";
  renderForensicPlot();

  $<HTMLButtonElement>("dlForensics").onclick = () => {
    downloadText(
      `forensic_parameters_${pop}.tsv`,
      tsvFromMatrix(headers, rows),
      "text/tab-separated-values",
    );
  };
  $<HTMLButtonElement>("dlForensicsXL").onclick = () => {
    writeXlsx(`forensic_parameters_${pop}.xlsx`, "forensics", headers, rows);
  };
}

function renderForensicPlot(): void {
  if (!lastForensicStats.length) return;
  const sel = $<HTMLSelectElement>("forensicCol");
  const key = sel.value as keyof LocusStats;
  const labels: string[] = [];
  const values: number[] = [];
  for (const s of lastForensicStats) {
    const v = s[key];
    if (typeof v === "number" && Number.isFinite(v)) {
      labels.push(s.locus);
      values.push(v);
    }
  }
  barPlot($("forensicPlot"), labels, values, sel.options[sel.selectedIndex]!.text);
}

// --- POPGEN tab ------------------------------------------------------------
$<HTMLInputElement>("displayDiv").addEventListener("change", (e) => {
  $("popgenWrap").style.display = (e.target as HTMLInputElement).checked ? "" : "none";
  if ((e.target as HTMLInputElement).checked) renderPopgen();
});
$<HTMLSelectElement>("popgenPop").addEventListener("change", renderPopgen);
$<HTMLSelectElement>("popgenCol").addEventListener("change", () => renderPopgenPlot());

interface PopgenRow {
  locus: string;
  N: number;
  Nall: number;
  Hobs: number | null;
  Hexp: number;
  HW_p: number | null;
  HW_chi: number | null;
  Ht: number | null;
  Fis: number | null;
  Fst: number | null;
}

let lastPopgenRows: PopgenRow[] = [];

function renderPopgen(): void {
  if (!state.genos) return;
  if (!$<HTMLInputElement>("displayDiv").checked) {
    $("popgenWrap").style.display = "none";
    return;
  }
  $("popgenWrap").style.display = "";

  const pop = $<HTMLSelectElement>("popgenPop").value || "all";
  const mask = maskFor(pop);
  const stats = locusIndices(state.genos, mask);
  const hwe =
    state.genos.ploidy === 2
      ? hweChiSquare(state.genos, mask)
      : null;

  // Per-locus F-statistics (only meaningful when more than one population
  // exists in the dataset and we're showing "all"). To match the original R
  // STRAF UI: Ht/Hs/Fis come from hierfstat::basic.stats; Fst comes from
  // hierfstat::wc (W&C θ̂). We mirror that split here.
  let fstats: LocusFStats[] | null = null;
  let bstats: LocusBasicStats[] | null = null;
  const allPops = uniquePopulations(state.genos);
  if (state.genos.ploidy === 2 && allPops.length > 1 && pop === "all") {
    try {
      fstats = perLocusFStats(state.genos);
    } catch {
      fstats = null;
    }
    try {
      bstats = basicStatsPerLocus(state.genos);
    } catch {
      bstats = null;
    }
  }

  const rows: PopgenRow[] = stats.map((s, i) => ({
    locus: s.locus,
    N: s.N,
    Nall: s.Nall,
    Hobs: s.Hobs ?? null,
    Hexp: s.GD,
    HW_p: hwe ? hwe[i]!.pValue : null,
    HW_chi: hwe ? hwe[i]!.chiSq : null,
    Ht: bstats ? bstats[i]!.Ht : null,
    Fis: bstats ? bstats[i]!.Fis : null,
    Fst: fstats ? fstats[i]!.Fst : null,
  }));
  lastPopgenRows = rows;

  type Col = { key: keyof PopgenRow; label: string; digits: number };
  const cols: Col[] = [
    { key: "N", label: "N", digits: 0 },
    { key: "Nall", label: "Nall", digits: 0 },
    { key: "Hobs", label: "Hobs", digits: 4 },
    { key: "Hexp", label: "Hexp", digits: 4 },
  ];
  if (state.genos.ploidy === 2) {
    cols.push(
      { key: "HW_chi", label: "HW χ²", digits: 3 },
      { key: "HW_p", label: "HW p", digits: 4 },
    );
  }
  if (fstats) {
    cols.push(
      { key: "Ht", label: "Ht", digits: 4 },
      { key: "Fis", label: "Fis", digits: 4 },
      { key: "Fst", label: "Fst", digits: 4 },
    );
  }
  const headers = ["Locus", ...cols.map((c) => c.label)];
  const tableRows: (string | number | null)[][] = rows.map((r) => [
    r.locus,
    ...cols.map((c) => {
      const v = r[c.key];
      return typeof v === "number" || v === null ? v : null;
    }),
  ]);
  renderTable($("popgenTable"), headers, tableRows, (v, col) => {
    if (col === 0) return v as string;
    if (v === null || v === undefined) return "";
    if (typeof v !== "number") return String(v);
    return fmtNum(v, cols[col - 1]!.digits);
  });

  const sel = $<HTMLSelectElement>("popgenCol");
  const prev = sel.value;
  sel.innerHTML = cols.map((c) => `<option value="${c.key}">${c.label}</option>`).join("");
  if (cols.some((c) => c.key === prev)) sel.value = prev;
  else sel.value = "Hobs";
  renderPopgenPlot();

  $<HTMLButtonElement>("dlPopgen").onclick = () => {
    downloadText(
      `popgen_${pop}.tsv`,
      tsvFromMatrix(headers, tableRows),
      "text/tab-separated-values",
    );
  };
  $<HTMLButtonElement>("dlPopgenXL").onclick = () => {
    writeXlsx(`popgen_${pop}.xlsx`, "popgen", headers, tableRows);
  };
}

function renderPopgenPlot(): void {
  if (!lastPopgenRows.length) return;
  const sel = $<HTMLSelectElement>("popgenCol");
  const key = sel.value as keyof PopgenRow;
  const labels: string[] = [];
  const values: number[] = [];
  for (const r of lastPopgenRows) {
    const v = r[key];
    if (typeof v === "number" && Number.isFinite(v)) {
      labels.push(r.locus);
      values.push(v);
    }
  }
  barPlot($("popgenPlot"), labels, values, sel.options[sel.selectedIndex]!.text);
}

// --- LD ---
$<HTMLButtonElement>("runLD").addEventListener("click", () => {
  state.ran.ld = true;
  $("ldWrap").style.display = "";
  $<HTMLButtonElement>("runLD").disabled = true;
  renderLd();
});

function renderLd(): void {
  if (!state.genos) return;
  if (state.genos.loci.length < 2) {
    $("ldTable").innerHTML = `<p class="muted">At least 2 loci are required.</p>`;
    return;
  }
  setStatus("Computing LD…");
  setComputing("ldWrap", "Computing pairwise LD…");
  deferCompute(() => {
    try {
      const results = pairwiseLd(state.genos!);
      const headers = ["Locus 1", "Locus 2", "N", "df", "χ²", "p-value"];
      const rows: (string | number | null)[][] = results.map((r) => [
        r.locus1,
        r.locus2,
        r.N,
        r.df,
        r.chiSq,
        r.pValue,
      ]);
      renderTable($("ldTable"), headers, rows, (v, col) => {
        if (col <= 1) return v as string;
        if (v === null || v === undefined || (typeof v === "number" && !Number.isFinite(v))) return "";
        if (typeof v !== "number") return String(v);
        if (col === 2 || col === 3) return Number.isInteger(v) ? String(v) : v.toFixed(0);
        return fmtNum(v, col === 5 ? 4 : 3);
      });

      // Heatmap of -log10(p) on the locus×locus square.
      const { loci, matrix } = ldMatrix(results, state.genos!.loci);
      renderLdHeatmap(loci, matrix);

      $<HTMLButtonElement>("dlLD").onclick = () => {
        downloadText("ld_pvalues.tsv", tsvFromMatrix(headers, rows), "text/tab-separated-values");
      };
      setStatus("LD test done.", "ok");
    } catch (err) {
      setStatus(`LD test failed: ${(err as Error).message}`, "error");
      $<HTMLButtonElement>("runLD").disabled = false;
      state.ran.ld = false;
    } finally {
      setComputing("ldWrap", null);
    }
  });
}

function renderLdHeatmap(loci: string[], matrix: (number | null)[][]): void {
  // Cell colour by −log10(p); cell label is the p-value itself.
  const values: (number | null)[][] = matrix.map((row) =>
    row.map((p) => (p === null || !Number.isFinite(p) ? null : -Math.log10(Math.max(p, 1e-10)))),
  );
  const labels: (string | null)[][] = matrix.map((row) =>
    row.map((p) => {
      if (p === null || !Number.isFinite(p)) return null;
      if (p < 1e-4) return "<.0001";
      return p.toFixed(3).replace(/^0/, "").replace(/^-0/, "-");
    }),
  );
  renderHeatmap($("ldHeatmap"), {
    values,
    labels,
    rowLabels: loci,
    colLabels: loci,
    vmin: 0,
    vmax: 5,
    legendTitle: "−log10 p",
    cellSize: 36,
  });
}

// --- Pairwise Fst ---
$<HTMLButtonElement>("runFst").addEventListener("click", () => {
  state.ran.fst = true;
  $("fstWrap").style.display = "";
  $<HTMLButtonElement>("runFst").disabled = true;
  renderFst();
});

function renderFst(): void {
  if (!state.genos) return;
  setStatus("Computing pairwise Fst…");
  setComputing("fstWrap", "Computing pairwise Fst…");
  deferCompute(() => {
    try {
      const fst = pairwiseFst(state.genos!);
      const headers = ["", ...fst.populations];
      const rows: (string | number | null)[][] = fst.populations.map((p, i) => [
        p,
        ...fst.populations.map((_, j) => (i === j ? null : i < j ? null : fst.matrix[i]![j]!)),
      ]);
      renderTable($("fstTable"), headers, rows, (v, col) => {
        if (col === 0) return v as string;
        if (v === null || v === undefined) return "";
        if (typeof v !== "number") return String(v);
        return fmtNum(v, 4);
      });

      const fullRows: (string | number | null)[][] = fst.populations.map((p, i) => [
        p,
        ...fst.matrix[i]!.map((v, j) => (i === j ? 0 : v)),
      ]);
      $<HTMLButtonElement>("dlFst").onclick = () => {
        downloadText("pairwise_fst.tsv", tsvFromMatrix(headers, fullRows), "text/tab-separated-values");
      };
      $<HTMLButtonElement>("dlFstXL").onclick = () => {
        writeXlsx("pairwise_fst.xlsx", "pairwise_fst", headers, fullRows);
      };
      setStatus("Pairwise Fst done.", "ok");
    } catch (err) {
      setStatus(`Pairwise Fst failed: ${(err as Error).message}`, "error");
      $<HTMLButtonElement>("runFst").disabled = false;
      state.ran.fst = false;
    } finally {
      setComputing("fstWrap", null);
    }
  });
}

// --- PCA + MDS tab ---------------------------------------------------------
$<HTMLButtonElement>("runPCA").addEventListener("click", () => {
  state.ran.pca = true;
  $("pcaWrap").style.display = "";
  $<HTMLButtonElement>("runPCA").disabled = true;
  computePca();
});

$<HTMLInputElement>("showEllipses").addEventListener("change", renderPcaPlot);
$<HTMLSelectElement>("pcX").addEventListener("change", renderPcaPlot);
$<HTMLSelectElement>("pcY").addEventListener("change", renderPcaPlot);
$<HTMLInputElement>("pcaScale").addEventListener("change", () => {
  if (state.pcaResult) computePca();
});

function computePca(): void {
  if (!state.genos) return;
  setStatus("Computing PCA…");
  setComputing("pcaWrap", "Computing PCA…");
  deferCompute(() => {
    try {
      const scale = $<HTMLInputElement>("pcaScale").checked;
      const res = pca(state.genos!, undefined, { scale });
      state.pcaResult = res;
      const nComp = Math.min(5, res.eigenvalues.length);
      const opts = Array.from({ length: nComp }, (_, i) => `<option value="${i}">${i + 1}</option>`).join("");
      $<HTMLSelectElement>("pcX").innerHTML = opts;
      $<HTMLSelectElement>("pcY").innerHTML = opts;
      if (nComp >= 2) $<HTMLSelectElement>("pcY").value = "1";
      renderPcaPlot();

      $<HTMLButtonElement>("dlPCA").onclick = () => {
        const hdr = ["ind", "pop", ...Array.from({ length: nComp }, (_, i) => `PC${i + 1}`)];
        const r: (string | number | null)[][] = res.individuals.map((ind, i) => [
          ind,
          res.populations[i]!,
          ...Array.from({ length: nComp }, (_, k) => res.scores[i]![k]!),
        ]);
        downloadText("pca_coordinates.tsv", tsvFromMatrix(hdr, r), "text/tab-separated-values");
      };

      $<HTMLButtonElement>("dlPCAEig").onclick = () => {
        const hdr = ["allele", ...Array.from({ length: nComp }, (_, i) => `PC${i + 1}`)];
        const r: (string | number | null)[][] = res.alleleLabels.map((lab, i) => [
          lab,
          ...Array.from({ length: nComp }, (_, k) => res.loadings[i]?.[k] ?? null),
        ]);
        downloadText(
          "pca_eigenvectors.tsv",
          tsvFromMatrix(hdr, r),
          "text/tab-separated-values",
        );
      };

      setStatus(
        `PCA done. PC1 explains ${(res.variance[0]! * 100).toFixed(1)}% of variance.`,
        "ok",
      );
    } catch (err) {
      setStatus(`PCA failed: ${(err as Error).message}`, "error");
      $<HTMLButtonElement>("runPCA").disabled = false;
      state.ran.pca = false;
    } finally {
      setComputing("pcaWrap", null);
    }
  });
}

function renderPcaPlot(): void {
  const res = state.pcaResult;
  if (!res) return;
  const xi = Number($<HTMLSelectElement>("pcX").value);
  const yi = Number($<HTMLSelectElement>("pcY").value);
  const showEllipses = $<HTMLInputElement>("showEllipses").checked;

  const groups = new Map<string, { name: string; x: number[]; y: number[]; text: string[] }>();
  for (let i = 0; i < res.scores.length; i++) {
    const pop = res.populations[i]!;
    let g = groups.get(pop);
    if (!g) {
      g = { name: pop, x: [], y: [], text: [] };
      groups.set(pop, g);
    }
    g.x.push(res.scores[i]![xi]!);
    g.y.push(res.scores[i]![yi]!);
    g.text.push(res.individuals[i]!);
  }
  const groupArr = Array.from(groups.values());

  let ellipses: { name: string; x: number[]; y: number[] }[] | undefined;
  if (showEllipses) {
    ellipses = [];
    for (const g of groupArr) {
      const e = confidenceEllipse(g.x, g.y, 0.95);
      if (e) ellipses.push({ name: g.name, x: e.x, y: e.y });
    }
  }

  const xPct = (res.variance[xi]! * 100).toFixed(1);
  const yPct = (res.variance[yi]! * 100).toFixed(1);
  scatterPlot(
    $("pcaPlot"),
    groupArr,
    `PC${xi + 1} (${xPct}%)`,
    `PC${yi + 1} (${yPct}%)`,
    ellipses,
  );
}

// MDS
$<HTMLButtonElement>("runMDS").addEventListener("click", () => {
  state.ran.mds = true;
  $("mdsWrap").style.display = "";
  $<HTMLButtonElement>("runMDS").disabled = true;
  computeMds();
});
$<HTMLSelectElement>("mdsDist").addEventListener("change", () => {
  if (cachedMds) computeMds();
});
$<HTMLSelectElement>("mdsX").addEventListener("change", renderMdsPlot);
$<HTMLSelectElement>("mdsY").addEventListener("change", renderMdsPlot);

interface CachedMds {
  populations: string[];
  points: number[][];
  variance: number[];
  distance: number[][];
  /** Kruskal stress (%) — only set when isoMDS was used. */
  stress?: number;
}
let cachedMds: CachedMds | null = null;

function computeMds(): void {
  if (!state.genos) return;
  const pops = uniquePopulations(state.genos);
  if (pops.length < 2) {
    setStatus("MDS requires at least 2 populations.", "error");
    $("mdsWrap").style.display = "none";
    cachedMds = null;
    state.ran.mds = false;
    $<HTMLButtonElement>("runMDS").disabled = false;
    return;
  }
  setStatus("Computing MDS…");
  setComputing("mdsWrap", "Computing MDS…");
  deferCompute(() => {
    try {
      const method = $<HTMLSelectElement>("mdsDist").value as DistanceMethod;
      const algo = $<HTMLSelectElement>("mdsAlgo").value;
      const dist = populationDistances(state.genos!, method);

      let points: number[][];
      let variance: number[];
      let stress: number | undefined;

      if (algo === "iso" && dist.populations.length >= 3) {
        const k = Math.min(5, dist.populations.length - 1);
        const iso = isoMds(dist.matrix, k);
        points = iso.points;
        // Per-axis "variance" — for non-metric MDS we report the proportion of
        // squared coordinate variance on each axis as a rough surrogate.
        const axisVar = new Array(k).fill(0) as number[];
        for (const row of points) for (let a = 0; a < k; a++) axisVar[a]! += (row[a] ?? 0) ** 2;
        const total = axisVar.reduce((s, v) => s + v, 0) || 1;
        variance = axisVar.map((v) => v / total);
        stress = iso.stress;
      } else {
        const mds = classicalMds(dist.matrix);
        points = mds.points;
        variance = mds.variance;
      }

      cachedMds = {
        populations: dist.populations,
        points,
        variance,
        distance: dist.matrix,
        stress,
      };

      const nAxes = Math.min(5, points[0]!.length);
      const opts = Array.from({ length: nAxes }, (_, i) => `<option value="${i}">${i + 1}</option>`).join("");
      $<HTMLSelectElement>("mdsX").innerHTML = opts;
      $<HTMLSelectElement>("mdsY").innerHTML = opts;
      if (nAxes >= 2) $<HTMLSelectElement>("mdsY").value = "1";
      renderMdsPlot();
      renderMdsTree();
      if (stress !== undefined) {
        $("mdsStress").textContent = `Kruskal stress: ${stress.toFixed(2)}%`;
      } else {
        $("mdsStress").textContent = "";
      }
      setStatus("MDS done.", "ok");
    } catch (err) {
      setStatus(`MDS failed: ${(err as Error).message}`, "error");
      $<HTMLButtonElement>("runMDS").disabled = false;
      state.ran.mds = false;
    } finally {
      setComputing("mdsWrap", null);
    }
  });
}

// Trigger recompute when algorithm changes too.
$<HTMLSelectElement>("mdsAlgo").addEventListener("change", () => {
  if (cachedMds) computeMds();
});

function renderMdsPlot(): void {
  const m = cachedMds;
  if (!m) return;
  const xi = Number($<HTMLSelectElement>("mdsX").value);
  const yi = Number($<HTMLSelectElement>("mdsY").value);
  const xs = m.points.map((p) => p[xi]!);
  const ys = m.points.map((p) => p[yi]!);
  const xPct = (m.variance[xi]! * 100).toFixed(1);
  const yPct = (m.variance[yi]! * 100).toFixed(1);
  labeledScatter(
    $("mdsPlot"),
    m.populations,
    xs,
    ys,
    `MDS axis ${xi + 1} (${xPct}%)`,
    `MDS axis ${yi + 1} (${yPct}%)`,
  );
}

function renderMdsTree(): void {
  const m = cachedMds;
  if (!m) return;
  if (m.populations.length < 2) return;
  const tree = hclust(m.distance, m.populations, "average");
  const seg = dendrogramSegments(tree);
  dendrogramPlot($("mdsTree"), seg.leafX, seg.leafLabels, seg.segments, seg.height);
}

// --- Reference MDS (STRidER) -----------------------------------------------
let refData: ReferenceFrequencies | null = null;

const refTabBtn = document.querySelector<HTMLButtonElement>('[data-subtab="refmds"]')!;
refTabBtn.addEventListener("click", async () => {
  if (!refData) await loadDefaultReference();
});

async function loadDefaultReference(): Promise<void> {
  setStatus("Loading reference frequencies…");
  try {
    const res = await fetch("./STRidER_frequencies_2024-09-24.csv");
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const text = await res.text();
    refData = parseReferenceCsv(text);
    populateRefPopList(refData.populations);
    setStatus(`Loaded ${refData.populations.length} reference populations.`, "ok");
    renderRefMds();
  } catch (err) {
    setStatus(`Could not load reference: ${(err as Error).message}`, "error");
  }
}

$<HTMLInputElement>("refFile").addEventListener("change", async (e) => {
  const file = (e.target as HTMLInputElement).files?.[0];
  if (!file) return;
  try {
    const text = await file.text();
    refData = parseReferenceCsv(text);
    populateRefPopList(refData.populations);
    setStatus(`Loaded custom reference: ${refData.populations.length} populations.`, "ok");
    renderRefMds();
  } catch (err) {
    setStatus(`Reference parse failed: ${(err as Error).message}`, "error");
  }
});

$<HTMLInputElement>("includeCurrentRef").addEventListener("change", renderRefMds);
$<HTMLSelectElement>("refPopsSelect").addEventListener("change", renderRefMds);
$<HTMLButtonElement>("refSelectAll").addEventListener("click", () => {
  const sel = $<HTMLSelectElement>("refPopsSelect");
  for (const o of Array.from(sel.options)) o.selected = true;
  renderRefMds();
});
$<HTMLButtonElement>("refSelectNone").addEventListener("click", () => {
  const sel = $<HTMLSelectElement>("refPopsSelect");
  for (const o of Array.from(sel.options)) o.selected = false;
  renderRefMds();
});

function populateRefPopList(pops: string[]): void {
  const sel = $<HTMLSelectElement>("refPopsSelect");
  sel.innerHTML = pops.map((p) => `<option value="${p}" selected>${p}</option>`).join("");
}

function selectedRefPops(): string[] {
  const sel = $<HTMLSelectElement>("refPopsSelect");
  return Array.from(sel.selectedOptions).map((o) => o.value);
}

function renderRefMds(): void {
  if (!refData) return;
  const selected = selectedRefPops();
  if (selected.length < 2) {
    $("refMdsPlot").innerHTML = `<p class="muted">Select at least 2 populations.</p>`;
    $("refCommonLoci").textContent = "";
    return;
  }
  // Filter to selected populations, drop columns missing in any selected pop.
  const filtered = dropMissingColumns(refData, selected);
  const includeCurrent = $<HTMLInputElement>("includeCurrentRef").checked;

  let labels: string[] = filtered.populations
    .map((p, i) => ({ p, i }))
    .filter(({ p }) => selected.includes(p))
    .map(({ p }) => p);
  let values: number[][] = filtered.populations
    .map((p, i) => ({ p, i }))
    .filter(({ p }) => selected.includes(p))
    .map(({ i }) => filtered.values[i]!.map((v) => (v === null ? 0 : v)));
  let commonLoci = uniqueLociFromColumns(filtered.columns);

  if (includeCurrent && state.genos) {
    // Compute one allele-frequency vector per population in the loaded dataset
    // (matching adegenet::genind2genpop behavior in the original R version).
    const perPop = currentDatasetVectorsPerPop(state.genos, filtered.columns);
    if (perPop.commonCols.length === 0) {
      $("refCommonLoci").textContent =
        "No alleles in common between the dataset and the reference — cannot include current data.";
    } else {
      // Restrict the existing reference rows to common columns first.
      const keepIdx: number[] = [];
      for (let c = 0; c < filtered.columns.length; c++) {
        if (perPop.commonCols.includes(filtered.columns[c]!)) keepIdx.push(c);
      }
      values = values.map((row) => keepIdx.map((c) => row[c]!));
      commonLoci = uniqueLociFromColumns(keepIdx.map((c) => filtered.columns[c]!));
      // Append one row per loaded population, prefixed so users can tell
      // their own pops apart from the reference set.
      for (const { pop, byCol } of perPop.byPop) {
        const vec = keepIdx.map((c) => byCol.get(filtered.columns[c]!) ?? 0);
        values.push(vec);
        labels.push(`★ ${pop}`);
      }
    }
  }

  $("refCommonLoci").textContent = `Loci used: ${commonLoci.join(", ")} (${commonLoci.length})`;

  // Nei distance + classical MDS.
  const D = neiDistanceFromMatrix(values);
  const mds = classicalMds(D);
  const xs = mds.points.map((p) => p[0] ?? 0);
  const ys = mds.points.map((p) => p[1] ?? 0);
  const xPct = ((mds.variance[0] ?? 0) * 100).toFixed(1);
  const yPct = ((mds.variance[1] ?? 0) * 100).toFixed(1);
  labeledScatter(
    $("refMdsPlot"),
    labels,
    xs,
    ys,
    `MDS axis 1 (${xPct}%)`,
    `MDS axis 2 (${yPct}%)`,
  );

  // Companion UPGMA tree on the same Nei distance matrix.
  if (labels.length >= 2) {
    const tree = hclust(D, labels, "average");
    const seg = dendrogramSegments(tree);
    dendrogramPlot($("refMdsTree"), seg.leafX, seg.leafLabels, seg.segments, seg.height);
  }
}

function uniqueLociFromColumns(cols: string[]): string[] {
  const set = new Set<string>();
  for (const c of cols) set.add(c.split("_")[0]!);
  return Array.from(set).sort();
}

/**
 * Project the currently-loaded genotype dataset into the reference's column
 * space, with one allele-frequency vector per population in the loaded data
 * (matching `adegenet::genind2genpop` in the original R version).
 *
 * Returns one entry per loaded population, plus the list of `<locus>_<allele>`
 * columns that exist in any of the per-pop tables AND in the reference.
 */
function currentDatasetVectorsPerPop(
  genos: Genotypes,
  refCols: string[],
): {
  byPop: { pop: string; byCol: Map<string, number> }[];
  commonCols: string[];
} {
  const pops = uniquePopulations(genos);
  const byPop: { pop: string; byCol: Map<string, number> }[] = [];
  const allPresent = new Set<string>();

  for (const pop of pops) {
    const mask = genos.populations.map((p) => p === pop);
    const byCol = new Map<string, number>();
    for (let l = 0; l < genos.loci.length; l++) {
      const counts = new Map<string, number>();
      let n = 0;
      for (let i = 0; i < genos.alleles[l]!.length; i++) {
        if (!mask[i]) continue;
        for (const a of genos.alleles[l]![i]!) {
          if (a === null) continue;
          counts.set(a, (counts.get(a) ?? 0) + 1);
          n++;
        }
      }
      if (n === 0) continue;
      for (const [a, c] of counts) {
        const col = `${genos.loci[l]}_${a}`;
        byCol.set(col, c / n);
        allPresent.add(col);
      }
    }
    byPop.push({ pop, byCol });
  }

  const commonCols = refCols.filter((c) => allPresent.has(c));
  return { byPop, commonCols };
}

// --- Haploid haplotype statistics -----------------------------------------
$<HTMLSelectElement>("haploPop").addEventListener("change", renderHaploStats);

function renderHaploSection(): void {
  const section = $("haploSection");
  if (!state.genos || state.genos.ploidy !== 1) {
    section.style.display = "none";
    return;
  }
  section.style.display = "";
  // Populate haploPop selector
  const sel = $<HTMLSelectElement>("haploPop");
  const pops = ["all", ...uniquePopulations(state.genos)];
  sel.innerHTML = pops.map((p) => `<option value="${p}">${p}</option>`).join("");
  renderHaploStats();
}

function renderHaploStats(): void {
  if (!state.genos || state.genos.ploidy !== 1) return;
  const pop = $<HTMLSelectElement>("haploPop").value || "all";
  let res;
  try {
    res = haplotypeStatsForPop(state.genos, pop);
  } catch (err) {
    setStatus(`Haplotype stats failed: ${(err as Error).message}`, "error");
    return;
  }

  // Haplotype table.
  const headers = ["Haplotype ID", "Haplotype", "Count", "Frequency"];
  const rows: (string | number | null)[][] = res.haplotypes.map((h) => [
    h.id,
    h.haplotype,
    h.count,
    h.frequency,
  ]);
  renderTable($("haploTable"), headers, rows, (v, col) => {
    if (col <= 1) return v as string;
    if (typeof v !== "number") return "";
    if (col === 2) return String(v);
    return fmtNum(v, 4);
  });

  $("haploDiv").textContent = `Haplotype diversity h = ${fmtNum(res.diversity, 4)}  (n = ${res.n}, distinct haplotypes = ${res.haplotypes.length})`;

  $<HTMLButtonElement>("dlHaplo").onclick = () => {
    downloadText(
      `haplotypes_${pop}.tsv`,
      tsvFromMatrix(headers, rows),
      "text/tab-separated-values",
    );
  };

  // Pairwise difference heatmap (only if ≥ 2 distinct haplotypes).
  if (res.haplotypes.length >= 2) {
    renderHaploDistHeatmap(
      res.haplotypes.map((h) => h.id),
      res.pairDiff,
    );
  } else {
    $("haploDistPlot").innerHTML = `<p class="muted small">At least 2 distinct haplotypes required for pairwise differences.</p>`;
  }
}

function renderHaploDistHeatmap(labels: string[], matrix: number[][]): void {
  renderHeatmap($("haploDistPlot"), {
    values: matrix,
    rowLabels: labels,
    colLabels: labels,
    legendTitle: "# differences",
    cellSize: 32,
    colorStops: [
      { t: 0, rgb: [255, 255, 255] },
      { t: 0.5, rgb: [253, 174, 107] },
      { t: 1, rgb: [230, 85, 13] },
    ],
  });
}

// (Graphical parameters now live in per-plot popovers — see ui/plotControls.ts)

// --- File conversion ------------------------------------------------------

function convertGuard(): Genotypes | null {
  if (!state.genos) {
    setStatus("Load a dataset first.", "error");
    return null;
  }
  return state.genos;
}

function dlConverted(name: string, content: string, mime = "text/plain"): void {
  downloadText(name, content, mime);
}

$<HTMLButtonElement>("dlGenepop").addEventListener("click", () => {
  const g = convertGuard();
  if (!g) return;
  try {
    dlConverted("straf2genepop.txt", toGenepop(g));
  } catch (err) {
    setStatus(`Genepop export failed: ${(err as Error).message}`, "error");
  }
});

$<HTMLButtonElement>("dlArlequin").addEventListener("click", () => {
  const g = convertGuard();
  if (!g) return;
  try {
    dlConverted("straf2arlequin.arp", toArlequin(g));
  } catch (err) {
    setStatus(`Arlequin export failed: ${(err as Error).message}`, "error");
  }
});

function selectedConvPop(): string {
  return $<HTMLSelectElement>("convPop").value || "all";
}

$<HTMLButtonElement>("dlFamilias").addEventListener("click", () => {
  const g = convertGuard();
  if (!g) return;
  const pop = selectedConvPop();
  try {
    dlConverted(`${pop}_Familias.csv`, toFamilias(g, pop));
  } catch (err) {
    setStatus(`Familias export failed: ${(err as Error).message}`, "error");
  }
});

$<HTMLButtonElement>("dlEuroformix").addEventListener("click", () => {
  const g = convertGuard();
  if (!g) return;
  const pop = selectedConvPop();
  try {
    dlConverted(`${pop}_Euroformix.csv`, toAlleleFreqCsv(g, "euroformix", pop), "text/csv");
  } catch (err) {
    setStatus(`Euroformix export failed: ${(err as Error).message}`, "error");
  }
});

$<HTMLButtonElement>("dlSTRmix").addEventListener("click", () => {
  const g = convertGuard();
  if (!g) return;
  const pop = selectedConvPop();
  try {
    dlConverted(`${pop}_STRmix.csv`, toAlleleFreqCsv(g, "strmix", pop), "text/csv");
  } catch (err) {
    setStatus(`STRmix export failed: ${(err as Error).message}`, "error");
  }
});

$<HTMLButtonElement>("dlLRmix").addEventListener("click", () => {
  const g = convertGuard();
  if (!g) return;
  const pop = selectedConvPop();
  try {
    dlConverted(`${pop}_LRmix.csv`, toAlleleFreqCsv(g, "lrmix", pop), "text/csv");
  } catch (err) {
    setStatus(`LRmix export failed: ${(err as Error).message}`, "error");
  }
});
