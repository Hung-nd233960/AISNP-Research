# SEA-JPT-CN Ancestry Inference Pipeline
# Run a single stage: make <target>
# Run full pipeline:  make all
# Remove intermediates: make clean

CONDA_ENV  := aisnp
PYTHON     := conda run -n $(CONDA_ENV) python3
NB_DIR     := notebooks/sea_jpt_cn

# Read genomes_data root from paths.yaml (or default)
GENOMES    := $(shell python3 -c \
  "import yaml; d=yaml.safe_load(open('paths.yaml')); print(d.get('genomes_data','data/1000genomes'))" \
  2>/dev/null || echo "data/1000genomes")

CACHE      := $(GENOMES)/cache
OUTPUTS    := $(GENOMES)/outputs

# ---------------------------------------------------------------------------
# Stamp files — touch these to mark a stage complete
# ---------------------------------------------------------------------------
STAMP_DIR       := .pipeline_stamps
$(STAMP_DIR):
	mkdir -p $@

STAMP_PREFILTER := $(STAMP_DIR)/prefilter.done   # 01 + 02 + 03
STAMP_FST       := $(STAMP_DIR)/04_fst.done
STAMP_STAT      := $(STAMP_DIR)/04a_stat.done
STAMP_ML_FST    := $(STAMP_DIR)/05b_ml_fst.done
STAMP_ML_BOTH   := $(STAMP_DIR)/05c_ml_both.done
STAMP_ML_STAT   := $(STAMP_DIR)/05a_ml_stat.done
STAMP_EVAL      := $(STAMP_DIR)/06_eval.done

# ---------------------------------------------------------------------------
.PHONY: all clean clean-stamps help

all: $(STAMP_EVAL)
	@echo "Pipeline complete."

help:
	@echo "Targets:"
	@echo "  all         Run the full pipeline end-to-end"
	@echo "  prefilter   01 hard filter → 02 situational filter → 03 VCF-to-matrix"
	@echo "  fst         04 FST selection + PCA"
	@echo "  stat        03b statistical SNP selection + 04b analysis"
	@echo "  ml          05a/b/c ML training (FST-only, FST+stat, stat-only)"
	@echo "  eval        06 model evaluation"
	@echo "  clean       Delete all pipeline intermediate files"
	@echo "  clean-stamps  Reset stamps only (re-run without deleting data)"

# ---------------------------------------------------------------------------
# Pre-filtering: 01 → 02 → 03
# 03 (VCF-to-matrix) is a shared dependency for both FST and statistical paths
# ---------------------------------------------------------------------------
prefilter: $(STAMP_PREFILTER)

$(STAMP_PREFILTER): $(STAMP_DIR)
	@echo "==> 01 Hard filtering..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/pre-filtering/01_hard_filtering.ipynb
	@echo "==> 02 Situational filtering..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/pre-filtering/02_situational_filtering.ipynb
	@echo "==> 03 VCF export + genotype matrix..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/pre-filtering/03_vcf_to_matrix.ipynb
	touch $@

# ---------------------------------------------------------------------------
# FST path: 04 FST selection + PCA
# ---------------------------------------------------------------------------
fst: $(STAMP_FST)

$(STAMP_FST): $(STAMP_PREFILTER)
	@echo "==> 04 FST selection and PCA..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/fst/04b_fst_and_pca.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Statistical path: 04a selection
# ---------------------------------------------------------------------------
stat: $(STAMP_STAT)

$(STAMP_STAT): $(STAMP_PREFILTER)
	@echo "==> 04a Statistical SNP selection..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/statistical_v1/04a_snp_selection.ipynb
	touch $@

# ---------------------------------------------------------------------------
# ML training: three variants
# ---------------------------------------------------------------------------
ml: $(STAMP_ML_FST) $(STAMP_ML_BOTH) $(STAMP_ML_STAT)

$(STAMP_ML_FST): $(STAMP_FST)
	@echo "==> 05b ML training (FST-only)..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/fst/05b_fst_only_training.ipynb
	touch $@

$(STAMP_ML_BOTH): $(STAMP_FST) $(STAMP_STAT)
	@echo "==> 05c ML training (FST + statistical)..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/statistical_v1/05c_fst_and_stat_training.ipynb
	touch $@

$(STAMP_ML_STAT): $(STAMP_STAT)
	@echo "==> 05a ML training (statistical-only)..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/statistical_v1/05a_stat_only_training.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Evaluation: 06
# ---------------------------------------------------------------------------
eval: $(STAMP_EVAL)

$(STAMP_EVAL): $(STAMP_ML_FST) $(STAMP_ML_BOTH) $(STAMP_ML_STAT)
	@echo "==> 06 Model evaluation..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/06_model_evaluation.ipynb
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/06b_stat_evaluation.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Housekeeping
# ---------------------------------------------------------------------------
clean-stamps:
	rm -rf $(STAMP_DIR)
	@echo "Stamps cleared. Next 'make all' will re-run all stages."

clean: clean-stamps
	rm -rf $(CACHE) $(OUTPUTS)
	@echo "cache/ and outputs/ deleted."

clean-cache: clean-stamps
	rm -rf $(CACHE)
	@echo "cache/ deleted (outputs preserved)."
