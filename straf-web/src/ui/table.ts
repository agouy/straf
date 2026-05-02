/** Render a 2D table into an element. Cells are formatted with `fmt`. */
export function renderTable<T>(
  container: HTMLElement,
  headers: string[],
  rows: T[][],
  fmt: (v: T, col: number) => string = (v) => String(v),
): void {
  const table = document.createElement("table");
  const thead = document.createElement("thead");
  const trh = document.createElement("tr");
  for (const h of headers) {
    const th = document.createElement("th");
    th.textContent = h;
    trh.appendChild(th);
  }
  thead.appendChild(trh);
  table.appendChild(thead);

  const tbody = document.createElement("tbody");
  for (const row of rows) {
    const tr = document.createElement("tr");
    for (let i = 0; i < row.length; i++) {
      const td = document.createElement("td");
      td.textContent = fmt(row[i] as T, i);
      tr.appendChild(td);
    }
    tbody.appendChild(tr);
  }
  table.appendChild(tbody);

  container.replaceChildren(table);
}

export function fmtNum(v: number | null | undefined, digits = 4): string {
  if (v === null || v === undefined || Number.isNaN(v)) return "";
  if (!Number.isFinite(v)) return v > 0 ? "∞" : "-∞";
  if (Number.isInteger(v)) return String(v);
  return v.toFixed(digits);
}

/** Trigger a download of arbitrary text data with a given filename. */
export function downloadText(filename: string, content: string, mime = "text/plain"): void {
  const blob = new Blob([content], { type: mime });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

export function tsvFromMatrix(
  headers: string[],
  rows: (string | number | null)[][],
): string {
  const lines = [headers.join("\t")];
  for (const row of rows) {
    lines.push(
      row
        .map((c) => (c === null || c === undefined ? "" : typeof c === "number" ? formatNumberForTsv(c) : c))
        .join("\t"),
    );
  }
  return lines.join("\n") + "\n";
}

function formatNumberForTsv(v: number): string {
  if (!Number.isFinite(v)) return "";
  if (Number.isInteger(v)) return String(v);
  return v.toString();
}
