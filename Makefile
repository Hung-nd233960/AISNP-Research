# AISNP Research Pipeline — Han / JPT / SEA (1000 Genomes Phase 3)
#
# make all          Run the full pipeline end-to-end
# make sweep        Re-run Stage 2/3 only (after clearing its cache)
# make results      Re-render all figures and tables (notebook 11)
# make clean-sweep  Delete Stage 2/3 cache so notebook 08 re-runs fully
# make clean        Delete all intermediate outputs
# make help         Show this message

CONDA_ENV  := aisnp
NB_EXEC    := conda run -n $(CONDA_ENV) jupyter nbconvert \
              --to notebook --execute --inplace \
              --ExecutePreprocessor.timeout=7200

NB_DIR     := notebooks/sea_jpt_cn

# Resolve genomes_data root from paths.yaml (or paths.local.yaml if present)
GENOMES    := $(shell python3 -c "\
import yaml, os; \
p = 'paths.local.yaml' if os.path.exists('paths.local.yaml') else 'paths.yaml'; \
print(yaml.safe_load(open(p)).get('genomes_data','data/1000genomes'))" \
  2>/dev/null || echo "data/1000genomes")

OUTPUTS    := $(GENOMES)/outputs
SWEEP_DIR  := $(OUTPUTS)/self_evaluation/08_unified_panel_sweep

# ---------------------------------------------------------------------------
# Stamp files — touch to mark a stage complete; delete to force re-run
# ---------------------------------------------------------------------------
STAMP_DIR        := .pipeline_stamps
$(STAMP_DIR):
	mkdir -p $@

S_PREFILTER  := $(STAMP_DIR)/prefilter.done    # 01 → 02 → 03
S_FST        := $(STAMP_DIR)/fst.done          # 04b + 05b
S_STAT       := $(STAMP_DIR)/stat.done         # 04a + 05a
S_FSTSTAT    := $(STAMP_DIR)/fst_stat.done     # 05c
S_SWEEP      := $(STAMP_DIR)/sweep.done        # 08 unified panel sweep
S_BENCHMARK  := $(STAMP_DIR)/benchmark.done    # 09 published panel comparison
S_OVERLAP    := $(STAMP_DIR)/overlap.done      # 10 panel overlap / rsID map
S_RESULTS    := $(STAMP_DIR)/results.done      # 11 figures + tables

# ---------------------------------------------------------------------------
.PHONY: all preprocess prefilter fst stat fst_stat sweep benchmark overlap results \
        nested compare-nested clean-nested clean-sweep clean-stamps clean help

all: $(S_RESULTS)
	@echo "Pipeline complete. Figures: $(OUTPUTS)/self_evaluation/11_results/"

help:
	@echo ""
	@echo "AISNP Pipeline — Han / JPT / SEA"
	@echo ""
	@echo "Targets:"
	@echo "  preprocess   Normalize + merge raw VCFs → PLINK2 pfile (run ONCE first)"
	@echo "  all          Full pipeline (01 → 11)"
	@echo "  prefilter    01 hard filter + 02 situational filter + 03 matrix"
	@echo "  fst          04b FST/PCA + 05b FST-only ML"
	@echo "  stat         04a stat selection + 05a stat ML"
	@echo "  fst_stat     05c FST+stat consensus ML"
	@echo "  sweep        08 unified 3-stage ML sweep (main)"
	@echo "  nested       08b leak-free nested CV (in-fold pool selection)"
	@echo "  compare-nested  nested CV vs frozen leaky baseline (table+figure)"
	@echo "  benchmark    09 published panel comparison"
	@echo "  overlap      10 panel rsID conversion + overlap"
	@echo "  results      11 figures + tables"
	@echo ""
	@echo "Housekeeping:"
	@echo "  clean-sweep  Delete Stage 2/3 cache (forces re-run of 08 from Stage 2)"
	@echo "  clean-stamps Reset all stamps (re-runs all notebooks, keeps outputs)"
	@echo "  clean        Delete stamps + all pipeline outputs"
	@echo ""

# ---------------------------------------------------------------------------
# Data preparation — run ONCE before anything else.
# Requires 1000 Genomes VCF files already downloaded.
# Usage: make preprocess VCF_DIR=/path/to/vcfs [THREADS=16] [MEMORY_MB=32000]
# ---------------------------------------------------------------------------
VCF_DIR    ?= $(GENOMES)/raw_vcf
VCF_PREFIX ?= allchr
THREADS    ?= 0
MEMORY_MB  ?=

preprocess:
	@echo "==> VCF preprocessing: normalize + convert + merge..."
	@echo "    VCF_DIR=$(VCF_DIR)  OUTPUT=$(VCF_PREFIX)  THREADS=$(THREADS)"
	bash scripts/vcf_preprocessing.sh "$(VCF_DIR)" "$(VCF_PREFIX)" "$(THREADS)" "$(MEMORY_MB)"
	@echo "Preprocessing complete. Output: $(GENOMES)/plink/$(VCF_PREFIX).{pgen,pvar,psam}"

# ---------------------------------------------------------------------------
# Stage: Pre-filtering (01 → 02 → 03)
# ---------------------------------------------------------------------------
prefilter: $(S_PREFILTER)

$(S_PREFILTER): | $(STAMP_DIR)
	@echo "==> 01 Hard filtering..."
	$(NB_EXEC) $(NB_DIR)/pre-filtering/01_hard_filtering.ipynb
	@echo "==> 02 Situational filtering (HWE, LD pruning)..."
	$(NB_EXEC) $(NB_DIR)/pre-filtering/02_situational_filtering.ipynb
	@echo "==> 03 VCF export + genotype matrix cache..."
	$(NB_EXEC) $(NB_DIR)/pre-filtering/03_vcf_to_matrix.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage: FST path (04b FST/PCA + 05b FST-only ML)
# ---------------------------------------------------------------------------
fst: $(S_FST)

$(S_FST): $(S_PREFILTER)
	@echo "==> 04b FST selection and PCA..."
	$(NB_EXEC) $(NB_DIR)/fst/04b_fst_and_pca.ipynb
	@echo "==> 05b FST-only ML training..."
	$(NB_EXEC) $(NB_DIR)/fst/05b_fst_only_training.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage: Statistical path (04a selection + 05a ML)
# ---------------------------------------------------------------------------
stat: $(S_STAT)

$(S_STAT): $(S_PREFILTER)
	@echo "==> 04a Statistical SNP selection (chi2, dAF, JSD)..."
	$(NB_EXEC) $(NB_DIR)/statistical/04a_snp_selection.ipynb
	@echo "==> 05a Stat-only ML training..."
	$(NB_EXEC) $(NB_DIR)/statistical/05a_stat_training.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage: FST+stat consensus (05c — needs both paths)
# ---------------------------------------------------------------------------
fst_stat: $(S_FSTSTAT)

$(S_FSTSTAT): $(S_FST) $(S_STAT)
	@echo "==> 05c FST+stat consensus ML training..."
	$(NB_EXEC) $(NB_DIR)/statistical/05c_fst_stat_training.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage: Unified 3-stage ML sweep (08)
# ---------------------------------------------------------------------------
sweep: $(S_SWEEP)

$(S_SWEEP): $(S_FST) $(S_STAT) $(S_FSTSTAT)
	@echo "==> 08 Unified panel sweep (3-stage, ~30 min)..."
	$(NB_EXEC) $(NB_DIR)/self_evaluation/08_unified_panel_sweep.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage: Published panel benchmarking (09)
# ---------------------------------------------------------------------------
benchmark: $(S_BENCHMARK)

$(S_BENCHMARK): $(S_SWEEP)
	@echo "==> 09 Published panel comparison..."
	$(NB_EXEC) $(NB_DIR)/self_evaluation/09_published_panel_comparison.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage: Panel rsID conversion + overlap (10)
# ---------------------------------------------------------------------------
overlap: $(S_OVERLAP)

$(S_OVERLAP): $(S_SWEEP)
	@echo "==> 10 Panel overlap and rsID map..."
	$(NB_EXEC) $(NB_DIR)/self_evaluation/10_panel_overlap.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage: Results — all figures and tables (11)
# ---------------------------------------------------------------------------
results: $(S_RESULTS)

$(S_RESULTS): $(S_BENCHMARK) $(S_OVERLAP)
	@echo "==> 11 Generating figures and tables..."
	$(NB_EXEC) $(NB_DIR)/self_evaluation/11_results.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage: Leak-free nested CV (08b) + comparison to the leaky baseline
#   Rebuilds the stat/FST/fst_stat pools INSIDE each CV fold (Option B: per-fold
#   plink2 --fst), so held-out samples never influence SNP nomination. See
#   scripts/nested_cv_sweep.py and reports/METHODOLOGY.md §2.3.
# ---------------------------------------------------------------------------
NESTED_DIR := $(OUTPUTS)/self_evaluation/08b_nested_cv_sweep
CONDA_PY   := conda run --no-capture-output -n $(CONDA_ENV) python -u

nested: $(NESTED_DIR)/nested_tierB_results.csv

$(NESTED_DIR)/nested_tierB_results.csv: scripts/nested_cv_sweep.py \
                                        scripts/nested_selection.py scripts/nested_fst.py
	@echo "==> 08b Leak-free nested CV sweep (Tier B full curve + Tier A committed N)..."
	$(CONDA_PY) scripts/nested_cv_sweep.py --tier both

compare-nested: nested
	@echo "==> Comparing leak-free nested CV vs frozen leaky baseline..."
	$(CONDA_PY) scripts/nested_comparison.py
	@echo "Comparison: $(NESTED_DIR)/comparison_*.{csv,png,md}"

# ---------------------------------------------------------------------------
# Housekeeping
# ---------------------------------------------------------------------------
clean-nested:
	rm -rf $(NESTED_DIR)
	@echo "Nested-CV outputs cleared. Run 'make compare-nested' to regenerate."

clean-sweep:
	rm -f $(SWEEP_DIR)/stage2_results.csv \
	      $(SWEEP_DIR)/stage2_agg.csv \
	      $(SWEEP_DIR)/stage2_best.csv \
	      $(SWEEP_DIR)/stage3_results.csv \
	      $(SWEEP_DIR)/our_panels_oofpreds.pkl
	rm -f $(S_SWEEP) $(S_BENCHMARK) $(S_RESULTS)
	@echo "Stage 2/3 cache cleared. Run 'make sweep' or 'make all' to re-run."

clean-stamps:
	rm -rf $(STAMP_DIR)
	@echo "All stamps cleared. Next 'make all' will re-execute all notebooks."

clean: clean-stamps
	rm -rf $(OUTPUTS)
	@echo "Stamps and outputs/ deleted."
