# Leaky vs leak-free nested CV — comparison

Accuracy is 5-fold CV (seed 42). Δ = nested − leaky (negative = leakage removed).

| N | config | leaky acc | nested Tier B | Δ | nested Tier A |
|---|--------|-----------|---------------|---|---------------|
| 5 | stat+LR+LR | 69.25% | 67.46% ± 2.50 | -1.78 | — |
| 10 | fst_stat+EN+SVM_Lin | 75.99% | 72.22% ± 2.66 | -3.76 | — |
| 15 | stat+EN+SVM_RBF | 80.55% | 78.37% ± 2.03 | -2.18 | — |
| 20 | stat+EN+LR | 84.52% | 83.53% ± 1.35 | -0.99 | — |
| 25 | stat+EN+RF | 86.50% | 84.12% ± 2.27 | -2.38 | — |
| 30 | stat+EN+LR | 89.30% | 86.10% ± 2.56 | -3.19 | — |
| 35 | stat+EN+LR | 92.26% | 86.10% ± 4.63 | -6.16 | 84.91% |
| 40 | stat+EN+LR | 91.66% | 86.70% ± 3.29 | -4.97 | — |
| 45 | stat+EN+LR | 92.85% | 87.89% ± 1.62 | -4.95 | — |
| 50 | stat+EN+SVM_RBF | 93.05% | 91.67% ± 1.20 | -1.39 | 86.70% |
| 55 | stat+EN+SVM_RBF | 95.43% | 92.66% ± 1.72 | -2.77 | — |
| 60 | stat+EN+SVM_RBF | 95.24% | 92.07% ± 2.07 | -3.17 | — |
| 65 | stat+LR+SVM_RBF | 93.65% | 89.88% ± 4.22 | -3.77 | — |
| 70 | stat+EN+SVM_RBF | 95.04% | 92.46% ± 2.03 | -2.57 | 90.67% |
| 75 | stat+EN+LR | 94.84% | 91.87% ± 1.31 | -2.97 | — |
| 80 | stat+EN+SVM_RBF | 94.64% | 93.06% ± 1.24 | -1.58 | — |
| 90 | fst_stat+EN+SVM_RBF | 96.03% | 92.66% ± 3.89 | -3.37 | — |
| 100 | fst_stat+EN+LR | 96.23% | 92.86% ± 1.69 | -3.37 | — |

**Mean Δ across N:** -3.07 points.
**Committed panels (N=35/50/70) mean Δ:** -3.38 points.

External panels carry no selection optimism, so the nested numbers are the fair comparison to Cai/Cao/Shi. See `comparison_accuracy_vs_n.png`.