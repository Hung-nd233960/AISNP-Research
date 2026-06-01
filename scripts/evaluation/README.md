# Part 2: Known AISNPs Comparison

This module compares statistically-selected SNPs (from Part 1) with known Ancestry-Informative SNPs (AISNPs) from research papers and commercial products.

## Workflow

```
rsID list (CSV/TXT)
       ↓
[06_rsid_to_bed.ipynb]
       ↓
BED files (chr, start, end, rsID)
       ↓
[07_bed_to_ml_matrix.ipynb]
       ↓
ML matrices (sample, pop, SNP genotypes)
       ↓
[08_known_aisnps_ml.ipynb]
       ↓
Performance comparison, confusion matrices, plots
```

## Directory Structure

```
scripts/part2/
├── __init__.py           # Module exports
├── rsid_utils.py         # rsID → BED conversion
├── bed_to_matrix.py      # BED → ML matrix extraction
├── ml_comparison.py      # ML model comparison utilities
└── README.md             # This file

scripts/notebooks/part2/
├── 06_rsid_to_bed.ipynb         # Convert rsID lists to BED
├── 07_bed_to_ml_matrix.ipynb    # Extract genotypes, create ML matrices
└── 08_known_aisnps_ml.ipynb     # ML comparison, confusion matrices, plots

data/known_aisnps/
└── *.csv, *.txt          # rsID lists from papers/products (add your files here)

output/part2/
├── *.bed                 # BED files
├── *_ml_matrix.csv       # ML-ready matrices
├── ml_comparison_results.csv
└── ml_comparison_report.txt

graphs/part2/
└── *.png                 # Comparison plots
```

## Adding New AISNP Sources

1. Add rsID file to `data/known_aisnps/`:
   - CSV format: One column with rsIDs (header optional)
   - TXT format: One rsID per line

2. Edit `06_rsid_to_bed.ipynb` → `SOURCES` list:

   ```python
   SOURCES = [
       {
           'name': 'your_source',
           'file': KNOWN_AISNPS_DIR / 'your_file.csv',
           'description': 'Description of the source'
       },
       ...
   ]
   ```

3. Run notebooks 06 → 07 → 08

## Helper Functions

### rsid_utils.py

- `load_rsid_list(file_path)` - Load rsIDs from CSV/TXT
- `rsid_to_bed(rsids, output_path)` - Convert rsIDs to BED format
- `batch_rsid_lookup(rsids)` - Query Ensembl API for coordinates

### bed_to_matrix.py

- `bed_to_snp_list(bed_path)` - Load SNP IDs from BED
- `extract_genotypes_from_pfile(pfile, snp_ids)` - Extract SNPs from PLINK
- `create_ml_matrix(pfile, snp_ids, samples_csv)` - Create ML-ready matrix

### ml_comparison.py

- `run_model_comparison(df, classifiers)` - Run K-fold CV on all models
- `generate_confusion_matrices(df, classifiers)` - Generate confusion matrices
- `plot_performance_comparison(results_df)` - Create comparison plots
- `create_summary_report(results_df)` - Generate text report

## Notes

- Uses GRCh37 (hg19) coordinates by default to match 1000 Genomes Phase 3
- API lookups are cached to avoid repeated queries
- Statistical SNPs from Part 1 are automatically included for comparison
