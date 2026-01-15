# Part 2: Known AISNPs Comparison
# Helper modules for comparing statistically-selected SNPs with known AISNPs
# from research papers and commercial products

from .rsid_utils import (
    load_rsid_list,
    rsid_to_bed,
    query_rsid_coordinates,
    batch_rsid_lookup,
)

from .bed_to_matrix import (
    bed_to_snp_list,
    extract_genotypes_from_vcf,
    extract_genotypes_from_pfile,
    pfile_to_genotype_matrix,
    create_ml_matrix,
)

from .ml_comparison import (
    run_model_comparison,
    generate_confusion_matrices,
    plot_performance_comparison,
)

__all__ = [
    "load_rsid_list",
    "rsid_to_bed",
    "query_rsid_coordinates",
    "batch_rsid_lookup",
    "bed_to_snp_list",
    "extract_genotypes_from_vcf",
    "extract_genotypes_from_pfile",
    "create_ml_matrix",
    "run_model_comparison",
    "generate_confusion_matrices",
    "plot_performance_comparison",
]
