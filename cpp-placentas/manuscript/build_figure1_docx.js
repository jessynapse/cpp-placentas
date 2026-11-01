// Builds an editable Figure 1 (participant flow) as a Word table: each box is
// a bordered cell, so text, widths and borders can be changed in Word.
// Counts come from outputs/flow_counts.csv (script 33).
// Run from the repository folder: node manuscript/build_figure1_docx.js

const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType,
  BorderStyle, AlignmentType, VerticalAlign,
} = require("docx");

const ROOT = path.resolve(__dirname, "..");
const rows = fs.readFileSync(path.join(ROOT, "outputs", "flow_counts.csv"), "utf8").trim().split("\n").slice(1);
const n = rows.map(r => Number(r.slice(r.lastIndexOf(",") + 1)));
const f = x => x.toLocaleString("en-US");

const FONT = "Calibri", SIZE = 24; // 12 pt (PPE figure labels)
const LINE = { style: BorderStyle.SINGLE, size: 8, color: "000000" };
const NONE = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const W = [4300, 800, 4260];

const p = (text, o = {}) => new Paragraph({
  alignment: o.align || AlignmentType.CENTER, spacing: { before: 0, after: 0 },
  indent: o.indent ? { left: 200 } : undefined,
  children: [new TextRun({ text, font: FONT, size: o.size || SIZE, bold: o.bold })],
});
const cell = (lines, w, boxed, o = {}) => new TableCell({
  width: { size: w, type: WidthType.DXA },
  verticalAlign: VerticalAlign.CENTER,
  borders: boxed ? { top: LINE, bottom: LINE, left: LINE, right: LINE } : { top: NONE, bottom: NONE, left: NONE, right: NONE },
  margins: { top: 80, bottom: 80, left: 120, right: 120 },
  children: lines.map(l => typeof l === "string" ? p(l, o) : l),
});
const boxRow = lines => new TableRow({ children: [cell(lines, W[0], true), cell([""], W[1], false), cell([""], W[2], false)] });
const exclRow = lines => new TableRow({ children: [
  cell(["↓"], W[0], false, { size: 28 }),
  cell(["→"], W[1], false, { size: 28 }),
  cell(lines.map((l, i) => p(l, { align: AlignmentType.LEFT, indent: i > 0 })), W[2], true),
] });

const table = new Table({
  width: { size: W.reduce((a, b) => a + b, 0), type: WidthType.DXA }, columnWidths: W,
  rows: [
    boxRow(["Offspring in the Collaborative Perinatal Project", `n = ${f(n[0])}`, "(58,760 pregnancies, 48,197 women)"]),
    exclRow([`Excluded: no placental pathology record (n = ${f(n[1])})`]),
    boxRow(["Offspring with a placental pathology record", `n = ${f(n[2])}`]),
    exclRow([`Excluded (n = ${f(n[3] + n[4] + n[5])})`, `Multiple gestation (n = ${f(n[3])})`, `Missing MVM score (n = ${f(n[4])})`, `Missing AI (n = ${f(n[5])})`]),
    boxRow(["Analytic cohort: singleton pregnancies with MVM and AI data", `n = ${f(n[6])}`]),
    exclRow([`IH outcome not observed (n = ${f(n[7] + n[8])})`, `Died before 12 months (n = ${f(n[7])})`, `No 12-month IH examination (n = ${f(n[8])})`, "Accounted for with inverse probability of censoring weights"]),
    boxRow(["Infants with an observed IH outcome at 12 months", `n = ${f(n[9])}`, `(${f(n[10])} definite IH cases)`]),
  ],
});

const doc = new Document({
  styles: { default: { document: { run: { font: FONT, size: SIZE } } } },
  sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } } },
    children: [
      table,
      new Paragraph({ spacing: { before: 240 }, children: [new TextRun({ text: "FIGURE 1 ", font: FONT, size: 24, bold: true }),
        new TextRun({ text: "Flow of participants in the Collaborative Perinatal Project. AI, acute inflammation. IH, infantile haemangioma. MVM, maternal vascular malperfusion.", font: FONT, size: 24 })] }),
    ] }],
});
Packer.toBuffer(doc).then(b => {
  fs.writeFileSync(path.join(ROOT, "manuscript", "Wong_PPE_IH_figure1_editable.docx"), b);
  console.log("Saved: manuscript/Wong_PPE_IH_figure1_editable.docx");
});
