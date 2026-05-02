/**
 * Minimal XLSX writer.
 *
 * Writes a single-sheet workbook from headers + a 2D array of cell values.
 * Numbers are emitted as numeric cells; strings as inline-string cells (no
 * shared-strings table — slightly larger files but trivially simpler code).
 *
 * The output is a real .xlsx (zipped Office Open XML) that opens cleanly in
 * Excel, LibreOffice, and Google Sheets. We rely on `fflate` for zipping;
 * everything else is handwritten string templates — total ~150 lines of code,
 * compared to ~700 KB for the SheetJS bundle.
 */

import { zipSync, strToU8 } from "fflate";
import { downloadText } from "./table.ts";

export type CellValue = string | number | null | undefined;

export function writeXlsx(
  filename: string,
  sheetName: string,
  headers: string[],
  rows: CellValue[][],
): void {
  const buf = buildXlsx(sheetName, headers, rows);
  // Trigger a binary download via Blob.
  // Copy into a fresh ArrayBuffer to satisfy strict BlobPart typing on TS 5.6+.
  const bytes = new Uint8Array(buf.byteLength);
  bytes.set(buf);
  const blob = new Blob([bytes.buffer], {
    type: "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
  });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
}

function buildXlsx(sheetName: string, headers: string[], rows: CellValue[][]): Uint8Array {
  const safeName = sheetName.replace(/[\\/?*:[\]]/g, "_").slice(0, 31) || "Sheet1";
  const allRows: CellValue[][] = [headers, ...rows];

  const sheetXml = renderSheet(allRows);
  const files: Record<string, Uint8Array> = {
    "[Content_Types].xml": strToU8(CONTENT_TYPES),
    "_rels/.rels": strToU8(ROOT_RELS),
    "xl/_rels/workbook.xml.rels": strToU8(WORKBOOK_RELS),
    "xl/workbook.xml": strToU8(renderWorkbook(safeName)),
    "xl/worksheets/sheet1.xml": strToU8(sheetXml),
  };
  return zipSync(files, { level: 6 });
}

function renderSheet(rows: CellValue[][]): string {
  const lines: string[] = [];
  lines.push('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>');
  lines.push(
    '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
  );
  lines.push("<sheetData>");
  for (let r = 0; r < rows.length; r++) {
    const row = rows[r]!;
    const rowNum = r + 1;
    lines.push(`<row r="${rowNum}">`);
    for (let c = 0; c < row.length; c++) {
      const v = row[c];
      if (v === null || v === undefined || v === "") continue;
      const ref = `${colLetter(c)}${rowNum}`;
      if (typeof v === "number" && Number.isFinite(v)) {
        lines.push(`<c r="${ref}"><v>${v}</v></c>`);
      } else {
        const s = escapeXml(String(v));
        lines.push(`<c r="${ref}" t="inlineStr"><is><t xml:space="preserve">${s}</t></is></c>`);
      }
    }
    lines.push("</row>");
  }
  lines.push("</sheetData>");
  lines.push("</worksheet>");
  return lines.join("");
}

function colLetter(col: number): string {
  let s = "";
  let n = col + 1;
  while (n > 0) {
    const rem = (n - 1) % 26;
    s = String.fromCharCode(65 + rem) + s;
    n = Math.floor((n - 1) / 26);
  }
  return s;
}

function escapeXml(s: string): string {
  return s
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;")
    .replace(/"/g, "&quot;")
    .replace(/'/g, "&apos;");
}

const CONTENT_TYPES = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
  <Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
</Types>`;

const ROOT_RELS = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>`;

const WORKBOOK_RELS = `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
</Relationships>`;

function renderWorkbook(sheetName: string): string {
  return `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
          xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <sheets><sheet name="${escapeXml(sheetName)}" sheetId="1" r:id="rId1"/></sheets>
</workbook>`;
}

// Re-export downloadText so call sites can import everything from one module.
export { downloadText };
