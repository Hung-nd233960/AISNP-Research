## Winner consistency (leak-free grid, N=5..100 x 5 folds)

Winners tallied over 90 (N, fold) cells; mean acc = mean over all folds & N.

**Pool** — fold-win rate (mean acc):
- stat: 62%  (84.3%)
- fst_stat: 13%  (82.0%)
- FST: 24%  (80.8%)

**Reductor** — fold-win rate (mean acc):
- EN: 62%  (84.8%)
- RF: 28%  (84.2%)
- LR: 10%  (78.1%)

**Classifier** — fold-win rate (mean acc):
- SVM_RBF: 37%  (83.6%)
- LR: 21%  (83.0%)
- RF: 16%  (82.7%)
- XGB: 11%  (82.1%)
- SVM_Lin: 8%  (81.8%)
- GBM: 8%  (81.0%)

**Top configurations** (by fold-win count):

| config | mean acc over N | top-1 folds | within-1SE of best (of N) |
|--------|-----------------|-------------|----------------------------|
| stat+EN+LR | 86.59% | 11/90 (12.2%) | 5/18 (27.8%) |
| FST+EN+SVM_RBF | 86.41% | 9/90 (10.0%) | 8/18 (44.4%) |
| stat+EN+SVM_RBF | 86.92% | 8/90 (8.9%) | 7/18 (38.9%) |
| stat+EN+XGB | 85.71% | 5/90 (5.6%) | 3/18 (16.7%) |
| stat+RF+SVM_RBF | 85.41% | 5/90 (5.6%) | 1/18 (5.6%) |
| FST+EN+LR | 85.21% | 5/90 (5.6%) | 3/18 (16.7%) |
| stat+LR+RF | 83.54% | 5/90 (5.6%) | 0/18 (0.0%) |
| stat+EN+RF | 86.14% | 4/90 (4.4%) | 3/18 (16.7%) |