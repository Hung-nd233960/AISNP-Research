"""
Backward-compatibility re-exports.

Notebooks and legacy callers that do `from utils import X` continue to work.
New code should import directly from plink / vcf / stats / io.
"""

from plink import run_plink2_command, count_variants, count_samples  # noqa: F401
from vcf import run_bcftools_command, genotype_to_numeric, vcf_to_numeric_matrix, add_population_labels  # noqa: F401
from stats import analyze_afreq, analyze_hardy, analyze_fst, merge_fst_top_variants  # noqa: F401
from fileio import ensure_dir, get_project_root, save_report, read_panel, extract_population_samples, variant_to_bed, variants_to_bed_file  # noqa: F401
