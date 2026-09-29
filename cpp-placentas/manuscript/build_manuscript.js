// Builds the manuscript, supplement, and response-to-coauthors Word files
// from manuscript/content.js, manuscript/responses.js, and outputs/*.csv.
// Run from the repository folder: node manuscript/build_manuscript.js
// Requires the docx npm package (npm install docx).

const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType,
  BorderStyle, AlignmentType, PageOrientation, ImageRun, PageBreak, Footer,
  PageNumber, LineNumberRestartFormat, HeadingLevel,
} = require("docx");

const ROOT = path.resolve(__dirname, "..");
const OUT = path.join(ROOT, "outputs");
const MS = path.join(ROOT, "manuscript");
const C = require("./content.js");
const R = require("./responses.js");

const FONT = "Times New Roman";
const BODY_SIZE = 24;   // 12 pt
const TABLE_SIZE = 16;  // 8 pt

// ---------------------------------------------------------------------------
// CSV
// ---------------------------------------------------------------------------
function readCSV(file) {
  const txt = fs.readFileSync(path.join(OUT, file), "utf8").replace(/\r/g, "");
  const rows = []; let row = [], cell = "", q = false;
  for (let i = 0; i < txt.length; i++) {
    const ch = txt[i];
    if (q) {
      if (ch === '"' && txt[i + 1] === '"') { cell += '"'; i++; }
      else if (ch === '"') q = false; else cell += ch;
    } else if (ch === '"') q = true;
    else if (ch === ",") { row.push(cell); cell = ""; }
    else if (ch === "\n") { row.push(cell); rows.push(row); row = []; cell = ""; }
    else cell += ch;
  }
  if (cell.length || row.length) { row.push(cell); rows.push(row); }
  return rows.filter(r => r.length > 1 || r[0] !== "");
}
const fmtNum = s => s.replace(/\b(\d{4,})\b/g, m => Number(m).toLocaleString("en-US"));
const ref1 = s => s.replace(/1\.00 \(ref\)/g, "1 REF").replace(/\(-/g, "(−").replace(/, -/g, ", −");
const lab = s => s.replace("MVM any vs none", "Any MVM").replace("AI any vs none", "Any AI").replace("MVM score, per 1 point", "MVM score, per point").replace("AI score, per 1 compartment", "AI score, per compartment");
const foetal = s => s.replace(/Fetal/g, "Foetal").replace(/fetal/g, "foetal");

// ---------------------------------------------------------------------------
// Citations: [@a, @b] -> [1,2] in order of first use
// ---------------------------------------------------------------------------
const order = [];
function cite(text) {
  return text.replace(/\[(@[^\]]+)\]/g, (_, inner) => {
    const keys = inner.split(",").map(k => k.trim().replace(/^@/, ""));
    const nums = keys.map(k => {
      if (!C.references[k]) throw new Error("Missing reference: " + k);
      if (!order.includes(k)) order.push(k);
      return order.indexOf(k) + 1;
    }).sort((a, b) => a - b);
    // compress runs: 1,2,3 -> 1-3
    const parts = []; let s = nums[0], p = nums[0];
    for (let i = 1; i <= nums.length; i++) {
      if (nums[i] === p + 1) { p = nums[i]; continue; }
      parts.push(s === p ? `${s}` : p === s + 1 ? `${s},${p}` : `${s}–${p}`);
      s = p = nums[i];
    }
    return "[" + parts.join(",") + "]";
  });
}

// ---------------------------------------------------------------------------
// Paragraph helpers
// ---------------------------------------------------------------------------
const run = (text, o = {}) => new TextRun({ text, font: FONT, size: o.size || BODY_SIZE, bold: o.bold, italics: o.italics, superScript: o.sup });
const para = (children, o = {}) => new Paragraph({
  children: Array.isArray(children) ? children : [run(children, o)],
  alignment: o.align, spacing: { line: o.line || 480, after: o.after ?? 0, before: o.before ?? 0 },
  indent: o.indent, keepNext: o.keepNext,
});
const heading = (text, level) => para([run(text, { bold: true, italics: level === 2 })], { before: 240, after: 0, keepNext: true });

// ---------------------------------------------------------------------------
// Tables: horizontal rules only (top, under header, bottom)
// ---------------------------------------------------------------------------
const NONE = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const LINE = { style: BorderStyle.SINGLE, size: 8, color: "000000" };
function table(header, rows, widths, opts = {}) {
  const total = widths.reduce((a, b) => a + b, 0);
  const cell = (text, w, o = {}) => {
    const lines = String(text).split("\n");
    const indent = /^ {2,}/.test(lines[0]);
    return new TableCell({
      width: { size: w, type: WidthType.DXA },
      borders: { top: o.top ? LINE : NONE, bottom: o.bottom ? LINE : NONE, left: NONE, right: NONE },
      margins: { top: 30, bottom: 30, left: 60, right: 60 },
      children: lines.map(l => new Paragraph({
        alignment: o.first ? AlignmentType.LEFT : AlignmentType.CENTER,
        indent: indent && o.first ? { left: 200 } : undefined,
        spacing: { line: 240 },
        children: [run(l.replace(/^ +/, ""), { size: TABLE_SIZE, bold: o.bold })],
      })),
    });
  };
  const trs = [new TableRow({
    tableHeader: true,
    children: header.map((h, j) => cell(h, widths[j], { top: true, bottom: true, bold: true, first: j === 0 })),
  })];
  rows.forEach((r, i) => {
    const last = i === rows.length - 1;
    const isSection = r.slice(1).every(x => x === "") && !/^ /.test(r[0]);
    trs.push(new TableRow({
      cantSplit: true,
      children: r.map((x, j) => cell(x, widths[j], { bottom: last, first: j === 0, bold: opts.boldSections && isSection && j === 0 })),
    }));
  });
  return new Table({ width: { size: total, type: WidthType.DXA }, columnWidths: widths, rows: trs });
}
const caption = (label, text) => para([run(label + " ", { bold: true, size: 20 }), run(text, { size: 20 })], { line: 276, after: 120, keepNext: true });
const note = text => para([run(text, { size: 16 })], { line: 240, before: 60 });

// split first data row ("", "n = ...") into the header
function headerWithN(rows) {
  const [hdr, nrow, ...rest] = rows;
  if (nrow[0] === "") return { header: hdr.map((h, j) => j === 0 ? h : `${h}\n${nrow[j]}`), rows: rest };
  return { header: hdr, rows: rows.slice(1) };
}

// ---------------------------------------------------------------------------
// MAIN TABLES
// ---------------------------------------------------------------------------
function table1() {
  const { header, rows } = headerWithN(readCSV("table1_corrected.csv"));
  const w = [2960, 1280, 1280, 1280, 1280, 1280];
  return [
    caption("Table 1.", "Characteristics of 45,268 singleton pregnancies in the Collaborative Perinatal Project, overall and by placental maternal vascular malperfusion (MVM) and acute inflammation (AI)."),
    table(header, rows, w, { boldSections: true }),
    note("Values are % (n) unless stated otherwise. Percentages exclude missing values, and the number missing is shown for each characteristic. Perinatal characteristics are shown for description only and were not included in adjusted models. BMI, body mass index. SD, standard deviation."),
  ];
}

function table2() {
  const t = readCSV("table2_main.csv");
  const get = lab => t.find(r => r[0] === lab).slice(2, 5);
  const REF = ["1 REF", "1 REF", "1 REF"];
  const rows = [
    ["Maternal vascular malperfusion", "", "", "", ""],
    ["  None", "1.6% (475/29,376)", ...REF],
    ["  Any", "2.0% (225/11,324)", ...get("Any MVM vs none")],
    ["  Score, per 1 point", "", ...get("MVM score, per 1 point")],
    ["Acute inflammation", "", "", "", ""],
    ["  None", "1.6% (519/31,820)", ...REF],
    ["  Any", "2.0% (181/8,880)", ...get("Any AI vs none")],
    ["  Score, per 1 compartment", "", ...get("AI score, per 1 compartment")],
    ["Joint MVM and AI", "", "", "", ""],
    ["  Neither", "1.6% (364/23,391)", ...REF],
    ["  MVM only", "1.8% (155/8,429)", ...get("MVM only")],
    ["  AI only", "1.9% (111/5,985)", ...get("AI only")],
    ["  Both MVM and AI", "2.4% (70/2,895)", ...get("Both MVM and AI")],
  ];
  return [
    caption("Table 2.", "Associations of placental maternal vascular malperfusion (MVM) and acute inflammation (AI) with definite infantile haemangioma (IH) by one year of age (700 cases among 40,700 infants)."),
    table(["Placental exposure", "IH, % (cases/N)", "Crude\nRR (95% CI)", "Primaryᵃ\nRR (95% CI)", "Secondaryᵇ\nRR (95% CI)"], rows, [2760, 1800, 1600, 1600, 1600], { boldSections: true }),
    note("Relative risks (RR) and 95% confidence intervals (CI) from modified Poisson regression with generalised estimating equations clustered on mother, pooled across 20 imputed datasets. ᵃAdjusted for maternal age, pre-pregnancy body mass index, race and ethnicity, education, income, marital status, smoking, parity, pre-pregnancy diabetes, chronic hypertension, infant sex, and study site, with stabilised inverse probability of censoring weights. ᵇAs for the primary model but without study site. MVM score range 0 to 8. AI score is the number of involved placental compartments (0 to 7)."),
  ];
}

// ---------------------------------------------------------------------------
// SUPPLEMENTARY TABLES
// ---------------------------------------------------------------------------
function supplement() {
  const out = [];
  const add = (label, cap, header, rows, widths, foot, opts) => {
    out.push(caption(label, cap), table(header, rows, widths, opts || {}));
    if (foot) out.push(note(foot));
    out.push(new Paragraph({ children: [new PageBreak()] }));
  };
  const MODEL_NOTE = "Primary: adjusted for 11 covariates and study site, with inverse probability of censoring weights. Secondary: as primary without study site. Definite IH, 700 cases among 40,700 infants, unless stated otherwise.";

  // S1
  let s = headerWithN(readCSV("tableS1_by_ih.csv"));
  add("Table S1.", "Placental pathology and maternal and infant characteristics among 40,700 infants with an observed outcome, by definite infantile haemangioma (IH) status.",
    s.header, s.rows.map(r => r.map(foetal)), [3960, 1800, 1800, 1800],
    "Values are % (n) unless stated otherwise. MVM, maternal vascular malperfusion. AI, acute inflammation. BMI, body mass index. SD, standard deviation.", { boldSections: true });

  // S2
  const s2 = readCSV("tableS2_site.csv");
  add("Table S2.", "Placental diagnoses, follow-up, and definite infantile haemangioma (IH) by Collaborative Perinatal Project study site (45,268 pregnancies).",
    s2[0], s2.slice(1), [2400, 1100, 1450, 1450, 1500, 1460],
    "Site names other than Boston are from the study team's site mapping and should be confirmed against CPP documentation. Cramér's V for the association of study site with any MVM, any AI, and definite IH among 40,700 infants: 0.38, 0.15, and 0.09. MVM, maternal vascular malperfusion. AI, acute inflammation.");

  // S3
  const s3 = readCSV("seq_cumulative_4exp.csv");
  const s3lab = { "1. Crude": "Crude", "2. + Demographics (age, race, education, income, marital)": "+ Demographic factorsᵃ",
    "3. + Health (BMI, smoking, diabetes, chronic HTN)": "+ Health factorsᵇ", "4. + Reproductive (parity, infant sex) = full, no site": "+ Parity and infant sex (all 11 covariates)",
    "5. + Study site = full with site": "+ Study site", "SECONDARY: full, no site, + IPW": "Secondary: 11 covariates + IPW", "PRIMARY: full + site + IPW": "Primary: 11 covariates + study site + IPW" };
  add("Table S3.", "Sequential adjustment of associations of placental pathology with definite infantile haemangioma.",
    ["Model", "Any MVM", "MVM score, per point", "Any AI", "AI score, per compartment"],
    s3.slice(1).map(r => [s3lab[r[0]] || r[0], ...r.slice(1)]), [3160, 1550, 1550, 1550, 1550],
    "Relative risks (95% CI). Covariate groups were added cumulatively without weights, followed by the secondary and primary models with inverse probability of censoring weights (IPW). ᵃMaternal age, race and ethnicity, education, income, marital status. ᵇPre-pregnancy body mass index, smoking, pre-pregnancy diabetes, chronic hypertension. MVM, maternal vascular malperfusion. AI, acute inflammation.");

  // S4
  const s4 = readCSV("mutual_4exp.csv");
  add("Table S4.", "Mutually adjusted associations of maternal vascular malperfusion (MVM) and acute inflammation (AI) with definite infantile haemangioma.",
    ["Exposure", "Mutually adjusted for", "Crude", "Primary", "Secondary", "Primary, single exposure"],
    s4.slice(1).map(r => [r[0], r[1], r[2], r[4], r[3], r[5]]), [2000, 1400, 1480, 1480, 1480, 1520], "Relative risks (95% CI). " + MODEL_NOTE);

  // S5
  const s5 = readCSV("synergy_reri.csv");
  add("Table S5.", "Additive and multiplicative interaction between maternal vascular malperfusion (MVM) and acute inflammation (AI) for definite infantile haemangioma.",
    ["Measure", "Crude", "Primary", "Secondary"], s5.slice(1).map(r => [r[0], r[1], r[2], r[3]].map(ref1)), [3960, 1800, 1800, 1800],
    "Estimates (95% CI) from the four-level joint exposure (reference: neither). CIs by the delta method using the pooled variance-covariance matrix. RERI, relative excess risk due to interaction. The synergy index is unstable when the single-exposure RRs are close to 1 and is shown for completeness. " + MODEL_NOTE);

  // S6
  const s6 = readCSV("grade_stage_separate.csv");
  add("Table S6.", "Associations of maternal vascular malperfusion (MVM) grade and acute inflammation (AI) stage with definite infantile haemangioma (IH).",
    ["Exposure", "IH, % (cases/N)", "Crude", "Primary", "Secondary"],
    s6.slice(1).flatMap((r, i, a) => {
      const row = [`  ${r[1]}`, r[2], r[5], r[6], r[7]].map(ref1);
      return (i === 0 || a[i - 1][0] !== r[0]) ? [[r[0], "", "", "", ""], row] : [row];
    }), [2360, 1800, 1700, 1700, 1800],
    "Relative risks (95% CI). MVM grade: low, score 1 to 2. High, score 3 or more. AI stage as recorded in the CPP placental dataset. " + MODEL_NOTE, { boldSections: true });

  // S7
  const s7 = readCSV("maternal_fetal_ai.csv"); const comp = readCSV("maternal_fetal_compartments.csv");
  add("Table S7.", "Associations of maternal and foetal acute inflammatory responses with definite infantile haemangioma (IH).",
    ["AI location", "IH, % (cases/N)", "Crude", "Primary", "Secondary"], s7.slice(1).map(r => r.slice(0, 5).map(foetal)), [2360, 1800, 1700, 1700, 1800],
    "Relative risks (95% CI). Maternal response: amnion or chorion of the membrane roll or placental surface. Foetal response: umbilical vein, umbilical artery, or foetal surface vessels. Among placentas with foetal-only inflammation (n = " + fmtNum(comp[2][1]) + "), the umbilical vein was involved in " + comp[2][6] + ", the umbilical artery in " + comp[2][7] + ", and foetal surface vessels in " + comp[2][8] + ". " + MODEL_NOTE);

  // S8
  const s8 = readCSV("boston_4exp.csv");
  add("Table S8.", "Associations of placental pathology with definite infantile haemangioma, stratified by the Boston centre and the other 11 centres.",
    ["Stratum", "Any MVM", "MVM score, per point", "Any AI", "AI score, per compartment"],
    s8.slice(1).map(r => [r[0].replace(/\[/, "\n(").replace(/\]/, ")").replace("Boston (site 5)", "Boston").replace("Other 11 sites (+ site)", "Other 11 centres").replace("All sites (primary)", "All centres (primary)")].concat(r.slice(1))),
    [2560, 1700, 1700, 1700, 1700],
    "Relative risks (95% CI), adjusted for 11 covariates with inverse probability of censoring weights, and for study site in the other 11 centres and all centres. In Boston, 92.6% of placentas were examined by one pathologist. MVM, maternal vascular malperfusion. AI, acute inflammation.");

  // S9
  const s9 = readCSV("boston_reader_validity.csv");
  add("Table S9.", "Associations of placental diagnoses with expected birth outcomes in the Boston centre and the other 11 centres, as a check of construct validity.",
    ["Outcome, by exposure", "Centre", "Exposed, % (n/N)", "Unexposed, % (n/N)", "Crude RR (95% CI)"],
    s9.slice(1).map((r, i) => [i % 2 === 0 ? r[0].replace("labor", "labour") : "", r[1].replace("sites", "centres"), fmtNum(r[2]), fmtNum(r[3]), r[4]]),
    [2560, 1500, 1800, 1800, 1700],
    "All 45,268 pregnancies with non-missing values. Pre-eclampsia includes superimposed pre-eclampsia. MVM, maternal vascular malperfusion. AI, acute inflammation. RR, relative risk.");

  // S10
  const s10 = readCSV("site_random_effect_models.csv");
  add("Table S10.", "Associations of placental pathology with definite infantile haemangioma, with study site as a fixed or random effect.",
    ["Exposure", "Crude", "Site as fixed effect (primary)", "Site as random effect", "No site (secondary)"],
    s10.slice(1).map(r => [lab(r[0]), r[3], r[4], r[5], r[6]]), [2560, 1700, 1700, 1700, 1700],
    "Relative risks (95% CI). Random-effect models were Poisson models with a random intercept for study site, adjusted for 11 covariates with inverse probability of censoring weights. Median rate ratio for study site: " + s10[1][7] + " to " + s10[4][7] + ". MVM, maternal vascular malperfusion. AI, acute inflammation.");

  // S11
  const s11 = readCSV("pathologist_models.csv");
  const ex11 = [...new Set(s11.slice(1).map(r => r[0]))];
  const pick = (e, m) => s11.find(r => r[0] === e && r[1] === m)[4];
  add("Table S11.", "Associations of placental pathology with definite infantile haemangioma, adjusted for the examining pathologist instead of study site.",
    ["Exposure", "Primary (study site)", "Examining pathologist", "Secondary (no site)"],
    ex11.map(e => [lab(e), pick(e, "Primary: covariates + site"), pick(e, "Covariates + examiner"), pick(e, "Secondary: covariates, no site")]),
    [3360, 2000, 2000, 2000],
    "Relative risks (95% CI), adjusted for 11 covariates with inverse probability of censoring weights. Examining pathologist was defined by site-specific examiner codes on the gross examination form, with examiners who read fewer than 100 placentas grouped within site. MVM, maternal vascular malperfusion. AI, acute inflammation.");

  // S12
  const s12 = readCSV("corrected_outcome_results.csv").filter(r => r[1] === "ih2");
  const ex12 = ["MVM any vs none", "MVM score (continuous)", "AI any vs none"];
  const p12 = (e, m) => (s12.find(r => r[2] === e && r[3] === m) || [])[7];
  add("Table S12.", "Associations of placental pathology with infantile haemangioma when suspect cases are classified as IH (930 cases among 40,700 infants).",
    ["Exposure", "Crude", "Primary", "Secondary"],
    ex12.map(e => [e.replace("MVM any vs none", "Any MVM").replace("MVM score (continuous)", "MVM score, per point").replace("AI any vs none", "Any AI"), p12(e, "crude"), p12(e, "adjusted+site+IPW"), p12(e, "adjusted+IPW")]),
    [3360, 2000, 2000, 2000], "Relative risks (95% CI). Prespecified sensitivity analysis. " + MODEL_NOTE.replace(" Definite IH, 700 cases among 40,700 infants, unless stated otherwise.", ""));

  // S13
  const f13 = readCSV("sibling_feasibility.csv"); const m13 = readCSV("sibling_models.csv");
  out.push(caption("Table S13.", "Feasibility and results of a within-mother (sibling) comparison."));
  out.push(table(["Step", "Mothers", "Pregnancies", "IH cases"], f13.slice(1), [4560, 1600, 1600, 1600]));
  out.push(para([run(" ", { size: 16 })], { line: 240 }));
  out.push(table(["Exposure", "Model", "OR (95% CI)"],
    m13.slice(1).map(r => [r[0], r[1].replace("Within-mother (conditional logistic)", "Within mother, conditional logistic").replace("Between-mother (logistic GEE)", "Between mothers, logistic GEE").replace(", same sibling sample", ""), r[2]]),
    [1800, 5560, 2000]));
  out.push(note("Among 11,900 pregnancies (181 IH cases) of 5,412 mothers with two or more infants in the sample. Only mothers discordant for both IH and the exposure contribute to within-mother estimates. Models adjusted for 11 covariates (and study site for between-mother models), pooled across 20 imputed datasets. OR, odds ratio. IH, infantile haemangioma."));
  return out;
}

// ---------------------------------------------------------------------------
// DOCUMENTS
// ---------------------------------------------------------------------------
const footer = new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 20 })] })] });
const page = { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } };
const brk = () => new Paragraph({ children: [new PageBreak()] });

function manuscript() {
  const kids = [];
  // title page
  kids.push(para([run(C.title, { bold: true, size: 28 })], { align: AlignmentType.CENTER, after: 240 }));
  const aut = [];
  C.authors.forEach(([n, s], i) => { aut.push(run(n)); aut.push(run(s, { sup: true })); if (i < C.authors.length - 1) aut.push(run(i === C.authors.length - 2 ? ", and " : ", ")); });
  kids.push(para(aut, { align: AlignmentType.CENTER, after: 240 }));
  C.affiliations.forEach(([s, t]) => kids.push(para([run(s, { sup: true, size: 20 }), run(t, { size: 20 })], { line: 276 })));
  kids.push(para([run("Corresponding author: ", { bold: true }), run(C.correspondence)], { before: 240 }));
  kids.push(para([run("Running head: ", { bold: true }), run(C.runningHead)]));
  kids.push(para([run("Manuscript type: ", { bold: true }), run("Original Article")]));
  const bodyWords = C.body.filter(b => b[0] === "p").map(b => b[1].replace(/\[@[^\]]+\]/g, "")).join(" ").split(/\s+/).length;
  const absWords = C.abstract.map(a => a[0] + " " + a[1]).join(" ").split(/\s+/).length;
  kids.push(para([run("Word count: ", { bold: true }), run(`${bodyWords.toLocaleString("en-US")} (main text), ${absWords} (abstract)`)]));
  kids.push(para([run("Tables: ", { bold: true }), run("2. Figures: 1. Supplementary material: 1 file (Tables S1 to S13). References: " + Object.keys(C.references).length)]));
  kids.push(para([run("ORCID iDs: ", { bold: true }), run(C.orcid.join(". "))]));
  kids.push(brk());

  // abstract, highlights, keywords
  kids.push(heading("Abstract", 1));
  C.abstract.forEach(([h, t]) => kids.push(para([run(h + ": ", { bold: true }), run(t)])));
  kids.push(para([run("Keywords: ", { bold: true }), run(C.keywords)], { before: 240 }));
  kids.push(heading("Highlights", 1));
  C.highlights.forEach(h => kids.push(para([run("• " + h)])));
  kids.push(brk());

  // body
  C.body.forEach(([type, text]) => {
    if (type === "h1") kids.push(heading(text, 1));
    else if (type === "h2") kids.push(heading(text, 2));
    else kids.push(para(cite(text), { indent: { firstLine: 480 } }));
  });

  // declarations
  C.declarations.forEach(([h, t]) => { kids.push(heading(h, 1)); kids.push(para(t)); });

  // references
  kids.push(heading("References", 1));
  order.forEach((k, i) => kids.push(para([run(`[${i + 1}] ${C.references[k]}`)], { line: 360, indent: { left: 480, hanging: 480 } })));
  const unused = Object.keys(C.references).filter(k => !order.includes(k));
  if (unused.length) console.warn("Unused references:", unused.join(", "));

  // figure legend
  kids.push(brk());
  kids.push(heading("Figure legend", 1));
  C.figureLegends.forEach(([l, t]) => kids.push(para([run(l + " ", { bold: true }), run(t)])));
  return { kids, bodyWords, absWords };
}

async function main() {
  const { kids, bodyWords, absWords } = manuscript();
  const tablesSection = [...table1(), brk(), ...table2(), brk(),
    caption("Fig. 1.", C.figureLegends[0][1]),
    new Paragraph({ children: [new ImageRun({ type: "png", data: fs.readFileSync(path.join(MS, "figure1_flow.png")), transformation: { width: 600, height: 554 } })] })];

  const styles = { default: { document: { run: { font: FONT, size: BODY_SIZE } } } };
  const doc = new Document({
    styles,
    sections: [
      { properties: { page, lineNumbers: { countBy: 1, restart: LineNumberRestartFormat.CONTINUOUS } }, footers: { default: footer }, children: kids },
      { properties: { page }, footers: { default: footer }, children: tablesSection },
    ],
  });
  fs.writeFileSync(path.join(MS, "Wong_Placenta_IH_manuscript.docx"), await Packer.toBuffer(doc));

  const supp = new Document({
    styles,
    sections: [{ properties: { page }, footers: { default: footer }, children: [
      para([run("Supplementary Material", { bold: true, size: 28 })], { align: AlignmentType.CENTER }),
      para([run(C.title, { italics: true })], { align: AlignmentType.CENTER, after: 240 }),
      brk(), ...supplement()] }],
  });
  fs.writeFileSync(path.join(MS, "Wong_Placenta_IH_supplement.docx"), await Packer.toBuffer(supp));

  const resp = new Document({
    styles,
    sections: [{ properties: { page: { ...page, size: { width: 12240, height: 15840, orientation: PageOrientation.LANDSCAPE } } }, footers: { default: footer }, children: [
      para([run("Response to co-author comments", { bold: true, size: 28 })], { after: 120 }),
      para([run(R.intro, { size: 20 })], { line: 276, after: 200 }),
      table(["#", "Co-author", "Comment", "Response", "Where"], R.rows.map((r, i) => [String(i + 1), ...r]), [400, 1300, 4300, 6000, 1960]),
      para([run(" ")]),
      para([run("Open questions for co-authors", { bold: true })], { line: 276, before: 200 }),
      ...R.questions.map(q => para([run("• " + q, { size: 20 })], { line: 276 })),
    ] }],
  });
  fs.writeFileSync(path.join(MS, "Response_to_coauthors.docx"), await Packer.toBuffer(resp));
  console.log(`Main text words: ${bodyWords}. Abstract words: ${absWords}. References: ${order.length}.`);
}
main();
