# SEA-JPT-CN Ancestry Inference Pipeline
# Run a single stage: make <target>
# Run full pipeline:  make all
# Remove intermediates: make clean

CONDA_ENV  := aisnp
PYTHON     := conda run -n $(CONDA_ENV) python3
PLINK      := conda run -n $(CONDA_ENV) plink2
NB_DIR     := scripts/notebooks/sea_jpt_cn

# Read genomes_data root from paths.yaml (or default)
GENOMES    := $(shell python3 -c \
  "import yaml; d=yaml.safe_load(open('paths.yaml')); print(d.get('genomes_data','data/1000genomes'))" \
  2>/dev/null || echo "data/1000genomes")

OUTPUT     := $(GENOMES)/output_sea_jpt_cn

# ---------------------------------------------------------------------------
# Stamp files — touch these to mark a stage complete
# ---------------------------------------------------------------------------
STAMP_DIR  := .pipeline_stamps
$(STAMP_DIR):
	mkdir -p $@

STAMP_PREFILTER   := $(STAMP_DIR)/01_prefilter.done
STAMP_SITUATIONAL := $(STAMP_DIR)/02_situational.done
STAMP_FST         := $(STAMP_DIR)/03_fst.done
STAMP_STAT        := $(STAMP_DIR)/03b_stat.done
STAMP_ML_FST      := $(STAMP_DIR)/04a_ml_fst.done
STAMP_ML_BOTH     := $(STAMP_DIR)/04b_ml_both.done
STAMP_ML_STAT     := $(STAMP_DIR)/04c_ml_stat.done
STAMP_EVAL        := $(STAMP_DIR)/05_eval.done

# ---------------------------------------------------------------------------
.PHONY: all clean clean-stamps help

all: $(STAMP_EVAL)
	@echo "Pipeline complete."

help:
	@echo "Targets:"
	@echo "  all              Run the full pipeline end-to-end"
	@echo "  prefilter        01+02  Hard filter + situational filter (pre-filtering notebooks)"
	@echo "  fst              03     FST calculation and SNP selection"
	@echo "  stat             03b    Statistical SNP analysis"
	@echo "  ml               04a/b/c Train ML models (FST-only, FST+stat, stat-only)"
	@echo "  eval             05     Model evaluation"
	@echo "  clean            Delete all pipeline intermediate files"
	@echo "  clean-stamps     Reset stage stamps only (force re-run without deleting data)"

# ---------------------------------------------------------------------------
# Stage 01+02: Pre-filtering (hard + situational filters)
# ---------------------------------------------------------------------------
prefilter: $(STAMP_PREFILTER)

$(STAMP_PREFILTER): $(STAMP_DIR)
	@echo "==> 01 Hard filtering..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/pre-filtering/01_hard_filtering.ipynb
	@echo "==> 02 Situational filtering..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/pre-filtering/02_situational_filtering.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage 03: FST selection + PCA
# ---------------------------------------------------------------------------
fst: $(STAMP_FST)

$(STAMP_FST): $(STAMP_PREFILTER)
	@echo "==> 03 FST selection and PCA..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/03_fst_and_pca.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage 03b: Statistical SNP selection
# ---------------------------------------------------------------------------
stat: $(STAMP_STAT)

$(STAMP_STAT): $(STAMP_PREFILTER)
	@echo "==> 02b Statistical SNP selection..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/02b_statistical_snp_selection.ipynb
	@echo "==> 03b Statistical SNP analysis..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/03b_statistical_snp_analysis.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage 04: ML training (three variants)
# ---------------------------------------------------------------------------
ml: $(STAMP_ML_FST) $(STAMP_ML_BOTH) $(STAMP_ML_STAT)

$(STAMP_ML_FST): $(STAMP_FST)
	@echo "==> 04a ML training (FST-only)..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/04a_fst_only_training.ipynb
	touch $@

$(STAMP_ML_BOTH): $(STAMP_FST) $(STAMP_STAT)
	@echo "==> 04b ML training (FST + statistical)..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/04b_fst_and_stat_training.ipynb
	touch $@

$(STAMP_ML_STAT): $(STAMP_STAT)
	@echo "==> 04c ML training (statistical-only)..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/04c_stat_only_training.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Stage 05: Evaluation
# ---------------------------------------------------------------------------
eval: $(STAMP_EVAL)

$(STAMP_EVAL): $(STAMP_ML_FST) $(STAMP_ML_BOTH) $(STAMP_ML_STAT)
	@echo "==> 05 Model evaluation..."
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/05_model_evaluation.ipynb
	$(PYTHON) -m jupyter nbconvert --to notebook --execute \
		--inplace $(NB_DIR)/05b_stat_evaluation.ipynb
	touch $@

# ---------------------------------------------------------------------------
# Housekeeping
# ---------------------------------------------------------------------------
clean-stamps:
	rm -rf $(STAMP_DIR)
	@echo "Stamps cleared. Next 'make all' will re-run all stages."

clean: clean-stamps
	rm -rf $(OUTPUT)
	@echo "Intermediates deleted."
