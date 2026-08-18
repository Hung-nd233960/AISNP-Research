# Documentation Plan — Report + Slides (leak-free version)

*Branch `nested-cv-leakfree`. Covers the two deliverables built from the leak-free
nested-CV rework: a long-form report and a short slide deck. This document is the
plan only — no report/slide content yet.*

**Format requirements (added):**
- **Report must be Word-compatible.** Written in Markdown, exported to `.docx` via
  Pandoc. This constrains figure format (PNG, not SVG — Pandoc's docx writer
  handles raster images reliably; SVG support is inconsistent) and forbids raw
  HTML in the source (doesn't survive conversion). Assumption: this requirement
  applies to the report, not the slide deck — Word isn't a slide medium, and the
  deck stays short-form Markdown. Flag if slides should also target a
  Word/PowerPoint-editable format.
- **Citations must be IEEE style.** Numbered in-text markers `[1]`, `[2]` in
  citation order of first appearance; reference list in IEEE reference format.
  Resolves the open citation-style question in §5 below.

---

## 0. Inputs already available (gathered before this plan)

| Source | What it gives us |
|---|---|
| `reports/SYSTEM_OVERVIEW.md` | Problem statement, data, pipeline architecture |
| `reports/METHODOLOGY.md` | Full methods, §2.3.4 leak-free nested CV, §2.3.1 leakage caveat |
| `reports/NESTED_CV_MASTER_REPORT.md` | Executive summary, leak/fix, honest performance, winner analysis, recommendation |
| `README.md` | Condensed key-results table (leaky vs leak-free), quick start |
| `results_archive/nested_leakfree_v1/08b_nested_cv_sweep/*.csv` | Raw numbers: `comparison_leaky_vs_nested.csv`, `winner_per_N.csv`, `winner_config_leaderboard.csv`, `nested_tierA_results.csv`, `nested_tierB_results.csv` |
| `results_archive/baseline_leaky_v1/` | Frozen leaky baseline (commit a803362) for before/after comparison |

**Gap identified:** all six committed figures in `graphs/` (`fig1_stage1_heatmap` …
`fig6_overlap_heatmap`) were generated 2026-06-20, **before** the leak fix — they
reflect the leaky pipeline. Only two leak-free figures exist so far
(`comparison_accuracy_vs_n.png`, `winner_consistency.png`, both in
`results_archive/nested_leakfree_v1/08b_nested_cv_sweep/`). §3 below plans which
figures need regenerating vs. which can be kept as "leaky, for contrast."

---

## 1. Two deliverables, two constraints

| | Report | Slides |
|---|---|---|
| Length | Liberal — can explain SNP/AISNP concepts, model background, full methodology, conclusions in depth | Short — headline numbers and figures only, no derivations |
| Audience | Reader who may not know what an AISNP panel or nested CV is | Audience with a few minutes; assumes a verbal walkthrough fills gaps |
| Format | Markdown source → **exported to `.docx` via Pandoc** (Word-compatible deliverable), paper-style structure | Markdown slide deck (Marp-style: `---` separated slides), one idea per slide — stays Markdown, no Word export |
| Citations | **IEEE style**: numbered `[1]`, `[2]` in-text, IEEE-format reference list, generated via Pandoc `--citeproc` + IEEE CSL + a maintained `.bib` file so numbering survives the docx export | Compact `[1]`-style markers only where a claim needs attribution; full list on a references slide |
| Primary source of truth | `NESTED_CV_MASTER_REPORT.md` + `METHODOLOGY.md`, expanded with background sections | Distilled from the report's §1 executive summary + §3/§4 tables |

Both deliverables are written **after** the plan below is confirmed, in the order:
content → writing → figures → insertion → citations. The same content plan feeds
both; slides subtract, they don't diverge.

---

## 2. Step 1 — Planning the content

### 2.1 Report content plan (long-form, section by section)

1. **Introduction / background** *(new writing — not in existing reports)*
   - What is a SNP; what makes a SNP "ancestry-informative" (AISNP)
   - Why within-EAS discrimination (Han/JPT/SEA) is hard — same super-population
   - Applications: forensics, pharmacogenomics, population genetics (already listed
     in `SYSTEM_OVERVIEW.md` §"Why it matters" — expand each with 2-3 sentences)
   - Prior work: Cai 2024, Cao 2022, Shi 2019 panels (cite properly — see §5)

2. **Data** — condense `SYSTEM_OVERVIEW.md` §Data + `METHODOLOGY.md` §2.1
   - 1000 Genomes Phase 3, 504 EAS samples, 3-group merge, quality filtering table

3. **Methods** — restructure `METHODOLOGY.md` in full, but reordered so the leak
   and its fix are presented as *the method*, not an afterthought:
   - 3.1 Candidate set construction (stat/FST/fst_stat pools) — explain each
     statistic in plain language (χ², JSD, AFD, Hudson F_ST) before the table
   - 3.2 The three-stage ML sweep (Stage 1 reductor/pool selection, Stage 2
     classifier eval, Stage 3 commitment) — explain what a "reductor" does and why
     ElasticNet/L1-LR/RF are the three tried
   - 3.3 **The leakage problem** — explain nested CV / leakage as a concept first
     (why fitting a feature selector on all labels before CV inflates estimates),
     *then* show where it happened in this pipeline (candidate pools built on all
     504 samples pre-CV)
   - 3.4 The fix — per-fold pool rebuilding (Option B), Tier A vs Tier B defined
     in plain language before the results table
   - 3.5 Models glossary (short, plain-language primer) — one paragraph each for
     Random Forest, XGBoost, Logistic Regression, SVM (RBF/linear), GBM,
     ElasticNet/L1-LR as reductors. This is the "explain some models" content the
     report has license to include; slides will drop this entirely.

4. **Results**
   - 4.1 Leaky vs leak-free accuracy, all N (table from `comparison_summary.md`,
     mean Δ −3.07pts, worst −6.16 at N=35)
   - 4.2 Honest performance at committed panels N=35/50/70 (Tier A/B/best-fixed)
   - 4.3 Winner analysis — reducer (EN dominant), pool (N-dependent stat→FST
     crossover ≈55–65), classifier (soft, SVM_RBF/LR close)
   - 4.4 Comparison to published panels (Cai-34, Cao-19, Shi-142) using leak-free
     numbers only

5. **Discussion / conclusions**
   - Reuse `NESTED_CV_MASTER_REPORT.md` §5 recommendation almost verbatim
   - Suggested narrative paragraph (parity + automation + transparency framing)
   - Limitations: committed SNP lists still selected on all 504 samples (only the
     accuracy *estimate* is corrected); config-selection optimism in "best fixed
     config" column; single-cohort (1000G EAS) generalizability

6. **Reproducibility** — `make nested`, `make compare-nested`, script list (from
   `NESTED_CV_MASTER_REPORT.md` §7), archive locations

7. **References** — see §5 below

### 2.2 Slides content plan (short-form, one deck)

Target ~12–15 slides. Each bullet below = roughly one slide.

1. Title
2. Problem statement (1 sentence: smallest SNP panel to tell Han/JPT/SEA apart)
3. Data (504 samples, 3 groups, 614,759 SNPs after filtering) — one table
4. Pipeline diagram (reuse `graphs/workflow_overview.svg` if still accurate, else
   regenerate — see §3.2)
5. **The leak** (one slide: candidate pools built on all samples before CV)
6. **The fix** (one slide: rebuild pools per fold, Tier A/B)
7. Headline number: leaky vs leak-free accuracy curve (figure)
8. Winner: ElasticNet reducer, N-dependent pool (figure: `winner_consistency.png`)
9. Committed panel results table (N=35/50/70, leaky vs leak-free)
10. Comparison to published panels (Cai/Cao/Shi) — leak-free only
11. Recommendation (pre-specify one classifier; report nested CV as primary)
12. Conclusion (parity + automation + transparency narrative, 1 sentence)
13. (Optional) Limitations, 3 bullets max
14. References (compact, slide-appropriate — first-author + year only, full list
    lives in the report)

No model glossary, no statistic derivations, no methodology sub-steps — those stay
in the report; slides state the conclusion each step produced.

---

## 3. Step 2 — Writing content

- Report: write section-by-section following §2.1 order. Each section pulls facts
  from the source docs in §0 — no new numbers are computed, only prose/explanation
  is added around existing verified results.
- Slides: write only after the report's Results/Discussion sections are finalized,
  since slide bullets are extracted from finished report prose, not drafted in
  parallel (avoids drift between the two documents).
- Flag every plain-language explanation (SNP background, model glossary, leakage
  concept) as report-only content when writing, so it's obvious what to cut when
  building the slide version.

## 4. Step 3 — Planning images and figures

**Word-compat constraint:** the report embeds **PNG only**, never SVG — Pandoc's
docx writer either rasterizes SVGs inconsistently or drops them depending on
version/engine. Every `.svg` figure in `graphs/` already has a `.png` sibling
(confirmed in the file listing), so this is a selection rule, not new rendering
work — except for any *new* figure (pool crossover plot, Tier A/B bar chart),
which must be exported as PNG from the start. Slides may keep using SVGs since
they don't go through the Word pipeline.

| Figure | Status | Action needed |
|---|---|---|
| `fig1_stage1_heatmap` | Leaky (Jun 20) | Regenerate from leak-free grid (`grid_results.csv`) or relabel explicitly as "Stage-1 sweep, leaky, shown for method illustration" |
| `fig2_accuracy_vs_n` | Leaky (Jun 20) | Superseded by `comparison_accuracy_vs_n.png` (already leak-free) — use that instead in both docs |
| `fig3_panel_comparison` | Leaky (Jun 20) | Regenerate against leak-free committed-panel numbers before use in Results §4.4 |
| `fig4_confusion` | Leaky (Jun 20) | Decide if confusion matrices are needed leak-free; if yes, regenerate; if this is illustrative only, caption as leaky |
| `fig5_classifier_robustness` | Leaky (Jun 20) | Cross-check against `winner_config_leaderboard.csv` classifier spread (~2.6pt); regenerate if numbers differ |
| `fig6_overlap_heatmap` | Not selection-optimism-affected (panel overlap is post-hoc) | Safe to reuse as-is |
| `workflow_overview.svg` | Architecture diagram, not data-dependent | Reuse; verify it already shows per-fold pool rebuilding (Option B) or annotate |
| `comparison_accuracy_vs_n.png` | Leak-free (new) | Use directly — headline figure for both report §4.1 and slide 7 |
| `winner_consistency.png` | Leak-free (new) | Use directly — report §4.3 and slide 8 |

**New figures likely needed** (not yet generated, check before writing):
- Pool crossover plot (stat vs FST vs fst_stat accuracy by N) — currently only in
  table form (`README.md` §3.2 of master report); a figure would strengthen §4.3
- Tier A vs Tier B vs best-fixed-config comparison at N=35/50/70 (bar chart) —
  currently table-only in master report §3.1

Action: before Step 4 (insertion), run the regeneration/creation pass so every
figure referenced by either document is confirmed leak-free-accurate or
explicitly labeled otherwise. This is real work (script runs), not just planning
— flag to user before executing since it may take nontrivial compute time.

## 5. Step 5 — Planning citations

Citation targets already named across existing docs, to be formalized into a
consistent reference list (BibTeX or plain numbered list, TBD when writing):

- 1000 Genomes Project Consortium (Phase 3 data source)
- PLINK2 (Chang et al., or the cog-genomics citation used in the field)
- Bhatia et al. 2013 — Hudson F_ST estimator (already named in `METHODOLOGY.md` §2.2.2)
- Cai et al. 2024 — EAS-specific 34-SNP panel
- Cao et al. 2022 — 19-SNP panel
- Shi et al. 2019 — 36/59/98/142-SNP nested panels
- scikit-learn (classifiers/reductors), XGBoost (Chen & Guestrin 2016) — standard
  tool citations for the model glossary section
- MyVariant.info — used for rsID mapping (panel overlap analysis)

**Citation style: decided — IEEE.** Numbered markers `[1]`, `[2]`, … assigned in
order of first appearance in the report; full reference list at the end in IEEE
format (`[n] Initials. Surname, "Title," Venue, vol. X, pp. Y–Z, Year.` — or the
appropriate IEEE dataset/software variant for non-paper sources like PLINK2 and
1000 Genomes).

Mechanics for the Word-compat pipeline:
- Maintain a single `reports/references.bib` (BibTeX) as source of truth for all
  entries listed above.
- Fetch the IEEE CSL style file (`ieee.csl`) once and keep it alongside the bib
  file in `reports/` for reproducible Pandoc builds.
- Report build: `pandoc report.md --citeproc --csl=ieee.csl --bibliography=references.bib -o report.docx`
  — citeproc numbers citations and generates the reference list automatically, so
  numbering stays correct even if citation order changes during editing.
- In-text citation syntax in the Markdown source uses Pandoc's `[@citekey]` form,
  not hand-typed `[1]` — hand-numbering would desync the moment a citation is
  reordered or a source is added.
- Slides: no citeproc pipeline (deck isn't exported to Word) — compact
  first-author+year markers stay as plain text, full list on the references
  slide, consistent in content with the report's IEEE list but not
  auto-generated.

## 6. Step 6 — Adding citations

- Report: insert `[@citekey]` markers at each claim site (methods sources,
  published panel comparisons, tool citations); citeproc generates the numbered
  IEEE in-text markers and reference list on build — do not hand-number.
- Slides: compact citation line under any slide referencing external panels
  (Cai/Cao/Shi) or third-party tools; full list on the final references slide.

---

## Execution order (confirmed)

1. Planning the content — this document (§1–§2)
2. Writing content — report first, full draft; then slides distilled from it
3. Planning images/figures — §3 table above, decide regenerate vs. reuse vs. new
4. Insert images and figures — after regeneration pass completes
5. Planning citations — §5 above (style decided: IEEE via Pandoc citeproc)
6. Adding citations — final pass on both documents
7. **Word-export verification** — build `report.docx` via the Pandoc command in
   §5 and open it to confirm: tables render as real Word tables (not code
   blocks), all figures embedded and captioned, citation numbering matches
   in-text order, headings map to Word heading styles for a navigable outline.
   Not optional — Markdown that looks correct can still convert badly (broken
   tables, missing images), so this is the acceptance check for the
   Word-compatibility requirement, not a formality.

Each step should be confirmed with the user before moving to the next, since
figure regeneration (step 3→4) involves running pipeline scripts, and step 7
requires Pandoc + the IEEE CSL file to be available in the environment.
**Verified:** Pandoc 3.7.0.2 is installed with docx output support — only the
`ieee.csl` file still needs fetching before step 7.
