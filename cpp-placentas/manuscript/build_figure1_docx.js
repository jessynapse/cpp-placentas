// Builds an editable Figure 1 (participant flow) as a Word table, laid out
// like Figure S1 of the AJE Advances CPP paper: main boxes on a vertical
// spine at left, exclusion and missing-data boxes at right joined to the
// spine by horizontal connectors, bold headings and bulleted counts.
// Every box is a bordered cell, so text, widths and borders stay editable.
// Counts: outputs/flow_counts.csv (script 33) and outputs/missingness.csv.
// Run from the repository folder: node manuscript/build_figure1_docx.js

const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType,
  BorderStyle, AlignmentType, VerticalAlign, HeightRule,
} = require("docx");

const ROOT = path.resolve(__dirname, "..");
const readCounts = file => fs.readFileSync(path.join(ROOT, "outputs", file), "utf8").trim().split("\n").slice(1);
const n = readCounts("flow_counts.csv").map(r => Number(r.slice(r.lastIndexOf(",") + 1)));
const miss = Object.fromEntries(readCounts("missingness.csv").map(r => { const [v, k] = r.split(","); return [v, Number(k)]; }));
const f = x => x.toLocaleString("en-US");
const pct = (k, d) => { const p = 100 * k / d; return k > 0 && p < 0.05 ? "<0.1%" : p.toFixed(1) + "%"; };

const FONT = "Calibri", SIZE = 24; // 12 pt (PPE figure labels)
const LINE = { style: BorderStyle.SINGLE, size: 8, color: "000000" };
const NONE = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const W = [1950, 1950, 450, 4850]; // spine left half, spine right half, connector, right boxes

const run = (text, o = {}) => new TextRun({ text, font: FONT, size: SIZE, bold: o.bold, italics: o.italics });
const nLine = (k, rest) => new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 0, after: 0 },
  children: [run("n", { italics: true }), run(` = ${f(k)}${rest ? " " + rest : ""}`)] });
const textLine = t => new Paragraph({ alignment: AlignmentType.CENTER, spacing: { before: 0, after: 0 }, children: [run(t)] });
const heading = t => new Paragraph({ spacing: { before: 0, after: 60 }, children: [run(t, { bold: true })] });
const bullet = (label, value) => new Paragraph({ spacing: { before: 0, after: 0 }, indent: { left: 400, hanging: 200 },
  children: [run("▪  "), run(label + ": ", { bold: true }), run(value)] });
const note = t => new Paragraph({ spacing: { before: 60, after: 0 }, children: [new TextRun({ text: t, font: FONT, size: SIZE, italics: true })] });

const borders = (o = {}) => ({ top: o.top ? LINE : NONE, bottom: o.bottom ? LINE : NONE, left: o.left ? LINE : NONE, right: o.right ? LINE : NONE });
const box = { top: true, bottom: true, left: true, right: true };
const cell = (children, w, b = {}, o = {}) => new TableCell({
  width: { size: w, type: WidthType.DXA }, columnSpan: o.span, rowSpan: o.rowSpan,
  verticalAlign: o.valign || VerticalAlign.CENTER, borders: borders(b),
  margins: { top: 100, bottom: 100, left: 140, right: 140 },
  children: children.length ? children : [new Paragraph("")],
});
const tr = (cells, h) => new TableRow({ cantSplit: true, height: h ? { value: h, rule: HeightRule.ATLEAST } : undefined, children: cells });

// Main box on the spine (spans the two spine columns)
const mainRow = lines => tr([cell(lines, W[0] + W[1], box, { span: 2 }), cell([], W[2]), cell([], W[3])]);
// Plain spine segment between boxes
const spineRow = h => tr([cell([], W[0], { right: true }), cell([], W[1]), cell([], W[2]), cell([], W[3])], h);
// Side box joined to the spine: two rows, connector on the bottom border of the first
const sideRows = (content, h) => [
  tr([cell([], W[0], { right: true }), cell([], W[1], { bottom: true }), cell([], W[2], { bottom: true }),
      cell(content, W[3], box, { rowSpan: 2, valign: VerticalAlign.TOP })], h),
  tr([cell([], W[0], { right: true }), cell([], W[1]), cell([], W[2])], h),
];

const nExcl = n[1] + n[3] + n[4] + n[5];
const elig = n[6];
const table = new Table({
  width: { size: W.reduce((a, b) => a + b, 0), type: WidthType.DXA }, columnWidths: W,
  rows: [
    mainRow([nLine(n[0], "offspring"), textLine(`${f(58760)} pregnancies`), textLine(`${f(48197)} women`)]),
    ...sideRows([
      heading(`Exclusion Criteria (n = ${f(nExcl)}, ${pct(nExcl, n[0])}):`),
      bullet("No placental pathology record", `n = ${f(n[1])} (${pct(n[1], n[0])})`),
      bullet("Multiple gestation", `n = ${f(n[3])} (${pct(n[3], n[0])})`),
      bullet("Missing MVM score", `n = ${f(n[4])} (${pct(n[4], n[0])})`),
      bullet("Missing AI", `n = ${f(n[5])} (${pct(n[5], n[0])})`),
    ], 600),
    mainRow([nLine(elig, "singleton pregnancies"), textLine("with MVM and AI data")]),
    ...sideRows([
      heading(`IH Outcome Not Observed Among ${f(elig)} Singleton Pregnancies (n = ${f(n[7] + n[8])}, ${pct(n[7] + n[8], elig)}):`),
      bullet("Died before 12 months", `n = ${f(n[7])} (${pct(n[7], elig)})`),
      bullet("No 12-month IH examination", `n = ${f(n[8])} (${pct(n[8], elig)})`),
      note("Accounted for with inverse probability of censoring weights"),
    ], 560),
    spineRow(200),
    ...sideRows([
      heading(`Missing Covariate Data Among ${f(elig)} Singleton Pregnancies:`),
      bullet("Pre-pregnancy BMI", `n = ${f(miss.bmi)} (${pct(miss.bmi, elig)})`),
      bullet("Income", `n = ${f(miss.income)} (${pct(miss.income, elig)})`),
      bullet("Education", `n = ${f(miss.educ)} (${pct(miss.educ, elig)})`),
      bullet("Parity", `n = ${f(miss.parity)} (${pct(miss.parity, elig)})`),
      bullet("Smoking", `n = ${f(miss.smoking)} (${pct(miss.smoking, elig)})`),
      bullet("Pre-pregnancy diabetes", `n = ${f(miss.dm)} (${pct(miss.dm, elig)})`),
      bullet("Chronic hypertension", `n = ${f(miss.chronic_htn)} (${pct(miss.chronic_htn, elig)})`),
      bullet("Marital status", `n = ${f(miss.marital)} (${pct(miss.marital, elig)})`),
      note("Imputed using multiple imputation (50 datasets)"),
    ], 1100),
    spineRow(200),
    mainRow([nLine(n[9], "infants examined"), textLine("for IH at 12 months"), textLine(`${f(n[10])} definite IH cases`)]),
  ],
});

const LEGEND = "Flow of participants in the Collaborative Perinatal Project. AI, acute inflammation. BMI, body mass index. IH, infantile haemangioma. MVM, maternal vascular malperfusion.";
const figureLegend = () => new Paragraph({ spacing: { before: 240 }, children: [
  new TextRun({ text: "FIGURE 1 ", font: FONT, size: 24, bold: true }), new TextRun({ text: LEGEND, font: FONT, size: 24 })] });

// Exported so build_ppe.js can place the same editable figure in the all-tables file.
module.exports = { figure1Table: () => table, figureLegend, LEGEND };
if (require.main === module) {
  const doc = new Document({
    styles: { default: { document: { run: { font: FONT, size: SIZE } } } },
    sections: [{ properties: { page: { size: { width: 12240, height: 15840 }, margin: { top: 1080, bottom: 1080, left: 1260, right: 1260 } } },
      children: [table, figureLegend()] }],
  });
  Packer.toBuffer(doc).then(b => {
    fs.writeFileSync(path.join(ROOT, "manuscript", "Wong_PPE_IH_figure1_editable.docx"), b);
    console.log("Saved: manuscript/Wong_PPE_IH_figure1_editable.docx");
  });
}
