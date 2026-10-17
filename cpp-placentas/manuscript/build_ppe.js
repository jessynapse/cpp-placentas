// Builds the Paediatric and Perinatal Epidemiology (PPE) manuscript and
// Supporting Information Word files from manuscript/content_ppe.js and the
// 50-imputation results in outputs/*_m50.csv.
// Run from the repository folder: node manuscript/build_ppe.js
// Requires the docx npm package (npm install docx).
//
// PPE format: Calibri 12 pt, double spaced, no indent, space between
// paragraphs, no line numbers, page numbers bottom centre, superscript
// citations after punctuation, Vancouver references, tables with horizontal
// rules only and no bold except headings, footnotes as superscript letters,
// reference categories as "1.00 (Reference)", supplementary tables as eTables.

const fs = require("fs");
const path = require("path");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell, WidthType,
  BorderStyle, AlignmentType, PageOrientation, ImageRun, PageBreak, Footer, PageNumber,
} = require("docx");

const ROOT = path.resolve(__dirname, "..");
const OUT = path.join(ROOT, "outputs");
const MS = path.join(ROOT, "manuscript");
const C = require("./content_ppe.js");

const FONT = "Calibri";
const BODY_SIZE = 24;   // 12 pt
const TABLE_SIZE = 18;  // 9 pt
const REFLAB = "1.00 (Reference)";

// ---------------------------------------------------------------------------
// CSV and label helpers
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
const fmtNum = s => String(s).replace(/\b(\d{4,})\b/g, m => Number(m).toLocaleString("en-US"));
const refFix = s => String(s).replace(/^1 REF$/, REFLAB).replace(/1\.00 \(ref\)/g, REFLAB)
  .replace(/\(-/g, "(−").replace(/, -/g, ", −").replace(/^-/, "−");
const lab = s => s.replace("MVM any vs none", "Any MVM").replace("AI any vs none", "Any AI")
  .replace("MVM score, per 1 point", "MVM score, per point").replace("AI score, per 1 compartment", "AI score, per compartment");
const sgaInfant = r => r.map(x => x === "Small for gestational age" ? "Small for gestational age infant" : x);
const foetal = s => s.replace(/Fetal/g, "Foetal").replace(/fetal/g, "foetal");
const siteName = s => s.replace(/^\d+ /, "");

// ---------------------------------------------------------------------------
// Text runs: superscript citations after punctuation, eTable numbering,
// and ^x markers for superscript footnote letters
// ---------------------------------------------------------------------------
const order = [];      // reference keys in order of first citation
const etOrder = [];    // eTable keys in order of first mention
const etNum = k => { if (!etOrder.includes(k)) etOrder.push(k); return etOrder.indexOf(k) + 1; };
const etText = s => s.replace(/\{eT:([a-z0-9]+)\}/g, (_, k) => "eTable " + etNum(k));

function citeNums(inner) {
  const nums = inner.split(",").map(k => k.trim().replace(/^@/, "")).map(k => {
    if (!C.references[k]) throw new Error("Missing reference: " + k);
    if (!order.includes(k)) order.push(k);
    return order.indexOf(k) + 1;
  }).sort((a, b) => a - b);
  const parts = []; let s = nums[0], p = nums[0];
  for (let i = 1; i <= nums.length; i++) {
    if (nums[i] === p + 1) { p = nums[i]; continue; }
    parts.push(s === p ? `${s}` : p === s + 1 ? `${s},${p}` : `${s}-${p}`);
    s = p = nums[i];
  }
  return parts.join(",");
}

const tr = (text, o = {}) => new TextRun({ text, font: FONT, size: o.size || BODY_SIZE, bold: o.bold, italics: o.italics, superScript: o.sup });

// "^a" -> superscript a
function supRuns(text, o = {}) {
  const out = [];
  String(text).split(/(\^[a-z])/).forEach(part => {
    if (/^\^[a-z]$/.test(part)) out.push(tr(part.slice(1), { ...o, sup: true }));
    else if (part) out.push(tr(part, o));
  });
  return out;
}

function bodyRuns(text, o = {}) {
  text = etText(text);
  const out = [];
  const re = /\s?\[(@[^\]]+)\]([.,:)]?)/g;
  let last = 0, m;
  while ((m = re.exec(text))) {
    if (m.index > last) out.push(tr(text.slice(last, m.index), o));
    if (m[2]) out.push(tr(m[2], o));
    out.push(tr(citeNums(m[1]), { ...o, sup: true }));
    last = re.lastIndex;
  }
  if (last < text.length) out.push(tr(text.slice(last), o));
  return out;
}

const para = (children, o = {}) => new Paragraph({
  children, alignment: o.align,
  spacing: { line: o.line || 480, after: o.after ?? 240, before: o.before ?? 0 },
  indent: o.indent, keepNext: o.keepNext,
});
const heading = (text, level) => para([tr(text, { bold: level < 3, italics: level === 3 })], { before: 240, after: 120, keepNext: true });
const brk = () => new Paragraph({ children: [new PageBreak()] });
const words = s => s.replace(/\[@[^\]]+\]/g, "").replace(/\{eT:[a-z0-9]+\}/g, "eTable").split(/\s+/).filter(Boolean).length;

// ---------------------------------------------------------------------------
// Tables: horizontal rules only (top, under header, bottom). Bold header only.
// ---------------------------------------------------------------------------
const NONE = { style: BorderStyle.NONE, size: 0, color: "FFFFFF" };
const LINE = { style: BorderStyle.SINGLE, size: 8, color: "000000" };
function table(header, rows, widths) {
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
        indent: indent && o.first ? { left: 220 } : undefined,
        spacing: { line: 240 },
        children: supRuns(l.replace(/^ +/, ""), { size: TABLE_SIZE, bold: o.bold }),
      })),
    });
  };
  const trs = [new TableRow({
    tableHeader: true,
    children: header.map((h, j) => cell(h, widths[j], { top: true, bottom: true, bold: true, first: j === 0 })),
  })];
  rows.forEach((r, i) => trs.push(new TableRow({
    cantSplit: true,
    children: r.map((x, j) => cell(x, widths[j], { bottom: i === rows.length - 1, first: j === 0 })),
  })));
  return new Table({ width: { size: total, type: WidthType.DXA }, columnWidths: widths, rows: trs });
}
const caption = (label, text) => para([tr(label + " ", { bold: true, size: 22 }), ...supRuns(text, { size: 22 })], { line: 276, after: 120, keepNext: true });
const note = text => para(supRuns(text, { size: 18 }), { line: 240, before: 60, after: 120 });

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
  return [
    caption("TABLE 1", "Characteristics of 45,268 singleton pregnancies in the Collaborative Perinatal Project, overall and by placental maternal vascular malperfusion (MVM) and acute inflammation (AI)."),
    table(header, rows.map(sgaInfant), [2960, 1280, 1280, 1280, 1280, 1280]),
    note("Values are % (n) unless stated otherwise. Percentages exclude missing values, and the number missing is shown for each characteristic. Perinatal characteristics are shown for description only and were not included in multivariable models. BMI, body mass index. SD, standard deviation."),
  ];
}

function table2() {
  const t = readCSV("table2_main_m50.csv");
  const get = l => t.find(r => r[0] === l).slice(2, 5).map(refFix);
  const REF = [REFLAB, REFLAB, REFLAB];
  const rows = [
    ["Joint MVM and AI", "", "", "", ""],
    ["  Neither", "1.6% (364/23,391)", ...REF],
    ["  MVM only", "1.8% (155/8,429)", ...get("MVM only")],
    ["  AI only", "1.9% (111/5,985)", ...get("AI only")],
    ["  Both MVM and AI", "2.4% (70/2,895)", ...get("Both MVM and AI")],
    ["Maternal vascular malperfusion", "", "", "", ""],
    ["  None", "1.6% (475/29,376)", ...REF],
    ["  Any", "2.0% (225/11,324)", ...get("Any MVM vs none")],
    ["  Score, per point^c", "", ...get("MVM score, per 1 point")],
    ["Acute inflammation", "", "", "", ""],
    ["  None", "1.6% (519/31,820)", ...REF],
    ["  Any", "2.0% (181/8,880)", ...get("Any AI vs none")],
    ["  Score, per compartment^d", "", ...get("AI score, per 1 compartment")],
  ];
  return [
    caption("TABLE 2", "Associations of placental maternal vascular malperfusion (MVM) and acute inflammation (AI) with definite infantile haemangioma (IH) by 12 months of age (700 cases among 40,700 infants)."),
    table(["Placental exposure", "IH, % (cases/N)", "Crude\nRR (95% CI)", "Primary^a\nRR (95% CI)", "Secondary^b\nRR (95% CI)"], rows, [2760, 1800, 1600, 1600, 1600]),
    note("Relative risks (RR) and 95% confidence intervals (CI) from modified Poisson regression with generalised estimating equations clustered on mother, pooled across 50 imputed datasets. ^aAdjusted for maternal age, pre-pregnancy body mass index, race and ethnicity, education, income, marital status, smoking, parity, pre-pregnancy diabetes, chronic hypertension, infant sex, and study site, with stabilised inverse probability of censoring weights. ^bAs for the primary model but without study site. ^cMVM score range 0 to 8. ^dNumber of involved placental compartments (0 to 7)."),
  ];
}

// ---------------------------------------------------------------------------
// eTABLES (numbered in order of first mention in the manuscript)
// ---------------------------------------------------------------------------
const MODEL_NOTE = "Primary: adjusted for 11 covariates and study site, with inverse probability of censoring weights. Secondary: as primary without study site. Pooled across 50 imputed datasets. Definite IH, 700 cases among 40,700 infants, unless stated otherwise. MVM, maternal vascular malperfusion. AI, acute inflammation. RR, relative risk. CI, confidence interval.";

const ETABLES = {
  incl: () => {
    const s = headerWithN(readCSV("etable_included_excluded.csv"));
    return { cap: "Characteristics of infants included in and excluded from the analysis of infantile haemangioma at 12 months (45,268 singleton pregnancies).",
      header: s.header, rows: s.rows.map(r => sgaInfant(r.map(x => siteName(x.replace(/^   (\d+ )/, "   ")))).map((x, j) => j === 0 && /^\s*\d+ /.test(x) ? x : x)),
      widths: [3000, 1700, 1700, 1700, 1700],
      foot: "Values are % (n) unless stated otherwise. Excluded infants either died before the 12-month examination or survived without an IH examination. MVM, maternal vascular malperfusion. AI, acute inflammation. BMI, body mass index. SD, standard deviation." };
  },
  byih: () => {
    const s = headerWithN(readCSV("tableS1_by_ih.csv"));
    return { cap: "Placental pathology and maternal and infant characteristics among 40,700 infants examined at 12 months, by definite infantile haemangioma (IH) status.",
      header: s.header, rows: s.rows.map(r => sgaInfant(r.map(foetal))), widths: [3960, 1800, 1800, 1800],
      foot: "Values are % (n) unless stated otherwise. MVM, maternal vascular malperfusion. AI, acute inflammation. BMI, body mass index. SD, standard deviation." };
  },
  site: () => {
    const raw = readCSV("tableS2_site.csv"); const rr = readCSV("site_table_percentages_m50.csv");
    const rows = raw.map((r, i) => i === 0 ? [...r, "Site RR (95% CI)^a"] :
      [...r, r[0] === "All sites" ? "" : refFix((rr.find(x => x[0] === r[0]) || [])[2] || "")]);
    return { cap: "Placental diagnoses, follow-up, and definite infantile haemangioma (IH) by Collaborative Perinatal Project study site (45,268 pregnancies).",
      header: rows[0], rows: rows.slice(1).map(r => [siteName(r[0]), ...r.slice(1)]), widths: [2000, 1000, 1250, 1250, 1300, 1300, 1260],
      foot: "^aRR of definite IH for each site compared with Boston, from the primary model for any MVM (adjusted for 11 covariates, with inverse probability of censoring weights). Site names other than Boston are from the study team's site mapping and should be confirmed against CPP documentation. MVM, maternal vascular malperfusion. AI, acute inflammation. RR, relative risk. CI, confidence interval." };
  },
  inter: () => {
    const s = readCSV("synergy_reri_m50.csv");
    return { cap: "Additive and multiplicative interaction between maternal vascular malperfusion (MVM) and acute inflammation (AI) for definite infantile haemangioma.",
      header: ["Measure", "Crude", "Primary", "Secondary"], rows: s.slice(1).map(r => [r[0], r[1], r[2], r[3]].map(refFix)), widths: [3960, 1800, 1800, 1800],
      foot: "Estimates (95% CI) from the four-category joint exposure (reference: neither). CIs by the delta method using the pooled variance-covariance matrix. RERI, relative excess risk due to interaction. The synergy index is unstable when the single-exposure RRs are close to 1 and is shown for completeness. " + MODEL_NOTE };
  },
  evalue: () => {
    const s = readCSV("evalues_m50.csv");
    return { cap: "E-values for associations of placental pathology with definite infantile haemangioma (primary model).",
      header: ["Exposure", "RR (95% CI)", "E-value, point estimate", "E-value, CI limit"],
      rows: s.slice(1).map(r => [r[0], r[1], Number(r[2]).toFixed(2), Number(r[3]).toFixed(2)]), widths: [3600, 2000, 1900, 1900],
      foot: "The E-value is the minimum strength of association, on the RR scale, that an unmeasured confounder would need to have with both the exposure and the outcome, beyond the measured covariates, to fully explain the observed association. The E-value for the CI is 1.00 when the CI includes 1. " + MODEL_NOTE };
  },
  mutual: () => {
    const s = readCSV("mutual_4exp_m50.csv");
    return { cap: "Mutually adjusted associations of maternal vascular malperfusion (MVM) and acute inflammation (AI) with definite infantile haemangioma.",
      header: ["Exposure", "Mutually adjusted for", "Crude", "Primary", "Secondary", "Primary, single exposure"],
      rows: s.slice(1).map(r => [lab(r[0]), r[1], r[2], r[4], r[3], r[5]]), widths: [2000, 1400, 1480, 1480, 1480, 1520],
      foot: "RR (95% CI). " + MODEL_NOTE };
  },
  grade: () => {
    const s = readCSV("grade_stage_separate_m50.csv");
    return { cap: "Associations of maternal vascular malperfusion (MVM) grade and acute inflammation (AI) stage with definite infantile haemangioma (IH).",
      header: ["Exposure", "IH, % (cases/N)", "Crude", "Primary", "Secondary"],
      rows: s.slice(1).flatMap((r, i, a) => {
        const row = [`  ${r[1]}`, r[2], r[5], r[6], r[7]].map(refFix);
        return (i === 0 || a[i - 1][0] !== r[0]) ? [[r[0], "", "", "", ""], row] : [row];
      }), widths: [2360, 1800, 1700, 1700, 1800],
      foot: "RR (95% CI). MVM grade: low, score 1 to 2. High, score 3 or more. AI stage as recorded in the CPP placental dataset. " + MODEL_NOTE };
  },
  mf: () => {
    const s = readCSV("maternal_fetal_ai_m50.csv"); const comp = readCSV("maternal_fetal_compartments.csv");
    return { cap: "Associations of acute inflammation (AI) by maternal or foetal involvement with definite infantile haemangioma (IH).",
      header: ["AI involvement", "IH, % (cases/N)", "Crude", "Primary", "Secondary"],
      rows: s.slice(1).map(r => r.slice(0, 5).map(foetal).map(refFix)), widths: [2360, 1800, 1700, 1700, 1800],
      foot: "RR (95% CI). Maternal response: amnion or chorion of the membrane roll or placental surface. Foetal response: umbilical vein, umbilical artery, or foetal surface vessels. Among placentas with foetal-only inflammation (n = " + fmtNum(comp[2][1]) + "), the umbilical vein was involved in " + comp[2][6] + ", the umbilical artery in " + comp[2][7] + ", and foetal surface vessels in " + comp[2][8] + ". " + MODEL_NOTE };
  },
  seq: () => {
    const s = readCSV("seq_cumulative_4exp_m50.csv");
    const L = { "1. Crude": "Crude", "2. + Demographics (age, race, education, income, marital)": "+ Demographic factors^a",
      "3. + Health (BMI, smoking, diabetes, chronic HTN)": "+ Health factors^b", "4. + Reproductive (parity, infant sex) = full, no site": "+ Parity and infant sex (all 11 covariates)",
      "5. + Study site = full with site": "+ Study site", "SECONDARY: full, no site, + IPW": "Secondary: 11 covariates + IPW", "PRIMARY: full + site + IPW": "Primary: 11 covariates + study site + IPW" };
    return { cap: "Sequential adjustment of associations of placental pathology with definite infantile haemangioma.",
      header: ["Model", "Any MVM", "MVM score, per point", "Any AI", "AI score, per compartment"],
      rows: s.slice(1).map(r => [L[r[0]] || r[0], ...r.slice(1)]), widths: [3160, 1550, 1550, 1550, 1550],
      foot: "RR (95% CI). Covariate groups were added cumulatively without weights, followed by the secondary and primary models with inverse probability of censoring weights (IPW). ^aMaternal age, race and ethnicity, education, income, marital status. ^bPre-pregnancy body mass index, smoking, pre-pregnancy diabetes, chronic hypertension. MVM, maternal vascular malperfusion. AI, acute inflammation." };
  },
  re: () => {
    const s = readCSV("site_random_effect_models_m50.csv");
    return { cap: "Associations of placental pathology with definite infantile haemangioma, with study site as a fixed or random effect.",
      header: ["Exposure", "Crude", "Site as fixed effect (primary)", "Site as random effect", "No site (secondary)"],
      rows: s.slice(1).map(r => [lab(r[0]), r[3], r[4], r[5], r[6]]), widths: [2560, 1700, 1700, 1700, 1700],
      foot: "RR (95% CI). Random-effect models were Poisson models with a random intercept for study site, adjusted for 11 covariates with inverse probability of censoring weights. Median rate ratio for study site: " + s[1][7] + " to " + s[4][7] + ". MVM, maternal vascular malperfusion. AI, acute inflammation." };
  },
  boston: () => {
    const s = readCSV("boston_4exp_m50.csv");
    return { cap: "Associations of placental pathology with definite infantile haemangioma, stratified by the Boston centre and the other 11 centres.",
      header: ["Stratum", "Any MVM", "MVM score, per point", "Any AI", "AI score, per compartment"],
      rows: s.slice(1).map(r => [r[0].replace(/\[/, "\n(").replace(/\]/, ")").replace("Boston (site 5)", "Boston").replace("Other 11 sites (+ site)", "Other 11 centres").replace("All sites (primary)", "All centres (primary)")].concat(r.slice(1))),
      widths: [2560, 1700, 1700, 1700, 1700],
      foot: "RR (95% CI), adjusted for 11 covariates with inverse probability of censoring weights, and for study site in the other 11 centres and all centres. In Boston, 92.6% of placentas were examined by one pathologist. MVM, maternal vascular malperfusion. AI, acute inflammation." };
  },
  ih2: () => {
    const s = readCSV("ih2_m50.csv");
    const ex = ["MVM any vs none", "MVM score (continuous)", "AI any vs none", "AI score (continuous)"].filter(e => s.some(r => r[2] === e));
    const p = (e, m) => (s.find(r => r[2] === e && r[3] === m) || [])[7] || "";
    return { cap: "Associations of placental pathology with infantile haemangioma when suspect cases are classified as IH (930 cases among 40,700 infants).",
      header: ["Exposure", "Crude", "Primary", "Secondary"],
      rows: ex.map(e => [e.replace("MVM any vs none", "Any MVM").replace("MVM score (continuous)", "MVM score, per point").replace("AI any vs none", "Any AI").replace("AI score (continuous)", "AI score, per compartment"), p(e, "crude"), p(e, "adjusted+site+IPW"), p(e, "adjusted+IPW")]),
      widths: [3360, 2000, 2000, 2000], foot: "RR (95% CI). Prespecified sensitivity analysis. " + MODEL_NOTE.replace(" Definite IH, 700 cases among 40,700 infants, unless stated otherwise.", "") };
  },
};

function supplement() {
  const unused = Object.keys(ETABLES).filter(k => !etOrder.includes(k));
  if (unused.length) console.warn("eTables not cited in the text:", unused.join(", "));
  const out = [];
  etOrder.forEach((k, i) => {
    if (!ETABLES[k]) throw new Error("Unknown eTable: " + k);
    const t = ETABLES[k]();
    out.push(caption(`eTable ${i + 1}.`, t.cap), table(t.header, t.rows, t.widths));
    if (t.foot) out.push(note(t.foot));
    if (i < etOrder.length - 1) out.push(brk());
  });
  return out;
}

// ---------------------------------------------------------------------------
// MANUSCRIPT
// ---------------------------------------------------------------------------
const footer = new Footer({ children: [new Paragraph({ alignment: AlignmentType.CENTER, children: [new TextRun({ children: [PageNumber.CURRENT], font: FONT, size: 20 })] })] });
const page = { size: { width: 12240, height: 15840 }, margin: { top: 1440, bottom: 1440, left: 1440, right: 1440 } };

function manuscript() {
  const kids = [];
  const bodyWords = C.body.filter(b => b[0] === "p").reduce((n, b) => n + words(b[1]), 0);
  const absWords = C.abstract.reduce((n, a) => n + words(a[0] + " " + a[1]), 0);
  const synWords = C.synopsis.reduce((n, a) => n + words(a[0] + " " + a[1]), 0);

  // title page
  kids.push(para([tr(C.title, { bold: true, size: 28 })], { align: AlignmentType.CENTER }));
  const aut = [];
  C.authors.forEach(([n, s], i) => { aut.push(tr(n)); aut.push(tr(s, { sup: true })); if (i < C.authors.length - 1) aut.push(tr(i === C.authors.length - 2 ? ", and " : ", ")); });
  kids.push(para(aut, { align: AlignmentType.CENTER }));
  C.affiliations.forEach(([s, t]) => kids.push(para([tr(s, { sup: true }), tr(t)], { after: 0 })));
  kids.push(para([tr("Correspondence: ", { bold: true }), tr(C.correspondence)], { before: 240 }));
  kids.push(para([tr("Running title: ", { bold: true }), tr(C.runningTitle)]));
  kids.push(para([tr("Word count: ", { bold: true }), tr(`${bodyWords.toLocaleString("en-US")} (main text), ${absWords} (abstract)`)]));
  kids.push(para([tr("Tables and figures: ", { bold: true }), tr("2 tables and 1 figure. Supporting Information: " + etOrderPreview() + " eTables.")]));
  kids.push(brk());

  // abstract
  kids.push(heading("ABSTRACT", 1));
  C.abstract.forEach(([h, t]) => kids.push(para([tr(h + ": ", { bold: true }), tr(t)])));
  kids.push(para([tr("Keywords: ", { bold: true }), tr(C.keywords)]));
  kids.push(para([tr("Word count: ", { bold: true }), tr(`${bodyWords.toLocaleString("en-US")} (main text)`)]));
  kids.push(brk());

  // synopsis and social media
  kids.push(heading("SYNOPSIS", 1));
  C.synopsis.forEach(([h, t]) => kids.push(para([tr(h + ": ", { bold: true }), tr(t)])));
  kids.push(heading("SOCIAL MEDIA QUOTE", 1));
  kids.push(para([tr(C.socialMedia)]));
  kids.push(brk());

  // body
  C.body.forEach(([type, text]) => {
    if (type === "h1") kids.push(heading(text, 1));
    else if (type === "h2") kids.push(heading(text, 2));
    else if (type === "h3") kids.push(heading(text, 3));
    else kids.push(para(bodyRuns(text)));
  });

  // declarations
  C.declarations.forEach(([h, t]) => { kids.push(heading(h, 1)); kids.push(para([tr(t)])); });

  // references
  kids.push(brk());
  kids.push(heading("REFERENCES", 1));
  order.forEach((k, i) => kids.push(para([tr(`${i + 1}. ${C.references[k]}`)], { indent: { left: 480, hanging: 480 } })));
  const unused = Object.keys(C.references).filter(k => !order.includes(k));
  if (unused.length) console.warn("Unused references:", unused.join(", "));

  return { kids, bodyWords, absWords, synWords };
}
// eTable count is only known after the body is processed, so count tokens up front
function etOrderPreview() {
  const keys = [];
  C.body.forEach(b => (b[1].match(/\{eT:([a-z0-9]+)\}/g) || []).forEach(t => { const k = t.slice(4, -1); if (!keys.includes(k)) keys.push(k); }));
  return keys.length;
}

async function main() {
  const { kids, bodyWords, absWords, synWords } = manuscript();
  const fig = fs.readFileSync(path.join(MS, "figure1_flow.png"));
  const tablesSection = [
    ...table1(), brk(), ...table2(), brk(),
    heading("FIGURE LEGEND", 1),
    ...C.figureLegends.map(([l, t]) => para([tr(l + " ", { bold: true }), tr(t)])),
    new Paragraph({ children: [new ImageRun({ type: "png", data: fig, transformation: { width: 600, height: 554 } })] }),
  ];
  const styles = { default: { document: { run: { font: FONT, size: BODY_SIZE } } } };

  const doc = new Document({ styles, sections: [
    { properties: { page }, footers: { default: footer }, children: kids },
    { properties: { page }, footers: { default: footer }, children: tablesSection },
  ] });
  fs.writeFileSync(path.join(MS, "Wong_PPE_IH_manuscript.docx"), await Packer.toBuffer(doc));

  const supp = new Document({ styles, sections: [{ properties: { page: { ...page, size: { width: 12240, height: 15840, orientation: PageOrientation.LANDSCAPE } } }, footers: { default: footer }, children: [
    para([tr("Supporting Information", { bold: true, size: 28 })], { align: AlignmentType.CENTER }),
    para([tr(C.title, { italics: true })], { align: AlignmentType.CENTER }),
    ...etOrder.map((k, i) => para([tr(`eTable ${i + 1}. `, { bold: true, size: 22 }), tr(ETABLES[k]().cap, { size: 22 })], { line: 276, after: 80 })),
    brk(), ...supplement()] }] });
  fs.writeFileSync(path.join(MS, "Wong_PPE_IH_supporting_information.docx"), await Packer.toBuffer(supp));

  const chars = C.socialMedia.length;
  console.log(`Main text: ${bodyWords} words (limit 3,500). Abstract: ${absWords} (limit 300). Synopsis: ${synWords} (limit 125). Social media: ${chars} characters (limit 280). Running title: ${C.runningTitle.length} characters (limit 40). References: ${order.length} (limit 60). eTables: ${etOrder.length}.`);
}
main();
