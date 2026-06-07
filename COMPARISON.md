# Published AISNP Panels — Comparison Targets

Goal: benchmark our selected SNP panels (FST / V2 stat, ± FST) against
**published, peer-reviewed AISNP panels** on the *same* CN/JPT/SEA task (504 samples,
1000 Genomes), using the same classifiers + 5-fold CV + same metrics
(accuracy / F1-weighted / MCC / ROC-AUC).

> **How to collect:** for each panel, drop a file of its SNPs into
> `data/published_panels/<key>.txt` (one rsID per line, or one b37 `chr:pos` per line).
> The comparison notebook will map them onto our genotype matrix
> (IDs look like `19:54792079[b37]G,T`), subset, and run the identical pipeline.
> Counts below are from the literature — **verify against the source** before quoting.

## Summary table

| key            | Panel                              | Reference (verify)                          | # SNPs | Scope                         | Relevance to CN/JPT/SEA | Collected? (file) |
|----------------|------------------------------------|---------------------------------------------|--------|-------------------------------|-------------------------|-------------------|
| `kidd_55`      | Kidd 55-AISNP                      | Kidd et al. 2014, FSI Genet.                | 55     | Global (7 regions)            | Coarse (continental)    | ☐                 |
| `seldin_128`   | Seldin/Kosoy 128-AISNP            | Kosoy et al. 2009, Hum Mutat                | 128    | Continental                   | Coarse                  | ☐                 |
| `pakstis_2019` | Pakstis expanded AISNP            | Pakstis et al. 2019, FSI Genet.             | ~115   | Global, finer                 | Some East-Asian signal  | ☐                 |
| `eurasiaplex`  | Eurasiaplex                        | Phillips et al. 2013, FSI Genet.            | 23     | Europe vs South Asia          | Low (not E-Asian focus) | ☐                 |
| `forenseq`     | ForenSeq DNA Signature (AISNP)    | Illumina MiSeq FGx (Kidd panel)             | 56     | Global (commercial)           | Coarse                  | ☐                 |
| `precision_id` | Precision ID Ancestry Panel       | Thermo Fisher (Seldin 123 + Kidd 55)        | 165    | Global (commercial)           | Coarse                  | ☐                 |
| `snpforid_34`  | SNPforID 34-plex                  | Phillips et al. 2007, FSI Genet.            | 34     | Identity (+ some ancestry)    | Low                     | ☐                 |
| `hsiao_lin_hwa`| East-Asian fine-scale panel       | **FILL: Hsiao/Lin/Hwa et al.**              | ?      | Within East Asia              | **High (target)**       | ☐

## Per-panel notes (fill SNP lists as you collect)

### kidd_55 — Kidd 55-AISNP

- The most-cited forensic AISNP set; basis of the ForenSeq AISNP block.
- Designed for *continental* resolution → expected to under-perform on the
  within-East-Asian CN/JPT/SEA split. Good as a "standard baseline".
- Source list: Kidd Lab / FROG-kb. File → `data/published_panels/kidd_55.txt`

### seldin_128 — Seldin/Kosoy 128-AISNP

- Continental ancestry; widely used. File → `seldin_128.txt`

### pakstis_2019

- Expanded universal AISNP set; finer than 55-panel. File → `pakstis_2019.txt`

### eurasiaplex

- Optimised for Europe vs South Asia; included for completeness, low E-Asian power.

### forenseq / precision_id

- Commercial kits. ForenSeq AISNP = Kidd 56; Precision ID = Seldin+Kidd merge (~165).
- Useful "what a forensic lab would actually use" baselines.

### hsiao_lin_hwa  ← **most relevant to our task**

- East-Asian / Taiwanese fine-scale AISNP study (verify exact citation + SNP list).
- This is the key head-to-head: a panel actually designed for within-East-Asia.

### cal_et_al

- Placeholder — add citation, scope, and SNP list.

## Comparison protocol (what the notebook will do)

1. Load each `data/published_panels/<key>.txt`.
2. Map panel SNPs → our matrix columns:
   - rsID lists: need rsID→b37 coord mapping (dbSNP) since our IDs are `chr:pos[b37]A,B`.
   - coord lists (`chr:pos`): match directly on position.
   - Report match rate (how many of the panel's SNPs exist in our LD-pruned matrix).
3. Subset the genotype matrix to the matched SNPs.
4. Run the **same** classifier suite + 5-fold stratified CV.
5. Report acc / F1 / MCC / ROC-AUC, side-by-side with our panels at the same N.
6. Fair-comparison caveats to state: (a) our panels are tuned on these very samples
   (in-sample optimism) while published panels are external → published panels are a
   *conservative* bar; (b) match rate < 100% means a published panel may be evaluated
   on fewer SNPs than its nominal size — always report the matched count.

## Open question to resolve

- rsID→coord mapping source (dbSNP build matching our b37 matrix). Needed for any
  panel published as rsIDs rather than coordinates.
