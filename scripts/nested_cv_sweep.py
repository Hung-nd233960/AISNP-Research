"""Leak-free nested cross-validation panel sweep (Option B).

Fixes the candidate-pool leakage in notebook 08: the statistical / FST / fst_stat
pools are rebuilt from *each outer training fold only* (via nested_selection +
nested_fst), so held-out samples never influence which SNPs are nominated.

Two tiers:
  Tier B  — isolate the leak. For each N use the *baseline-committed* config
            (panel, reductor, classifier) and re-estimate it with in-fold pools.
            Directly comparable to the frozen leaky Stage-2 numbers; the drop is
            the leakage.
  Tier A  — defensible headline. Fully nested: pools built on the outer-train,
            the config (panel, reductor, classifier) chosen by an INNER CV on the
            outer-train, then scored once on the untouched outer-test fold.

Outputs CSVs under outputs/self_evaluation/08b_nested_cv_sweep/.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (accuracy_score, f1_score, matthews_corrcoef,
                             roc_auc_score)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.svm import SVC
from xgboost import XGBClassifier

sys.path.insert(0, str(Path(__file__).resolve().parent))
from notebook_init import setup
from nested_fst import fst_pool_for_samples
from nested_selection import stat_pool_indices, fst_stat_pool_indices

RANDOM_STATE = 42
N_CV = 5
PANEL_SIZES = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100]
REDUCTORS = ["RF", "LR", "EN"]
PANELS = ["stat", "FST", "fst_stat"]


# --------------------------------------------------------------------------
# Reductor + classifier definitions — identical to notebook 08
# --------------------------------------------------------------------------

def fit_reductor(rname, X_tr, y_tr):
    """Return (per-feature importance, fitted StandardScaler). Mirrors nb 08."""
    scaler = StandardScaler()
    X_sc = scaler.fit_transform(X_tr)
    if rname == "RF":
        red = RandomForestClassifier(n_estimators=200, random_state=RANDOM_STATE, n_jobs=-1)
        red.fit(X_tr, y_tr)
        return red.feature_importances_, scaler
    if rname == "LR":
        red = LogisticRegression(max_iter=2000, solver="saga", penalty="l2",
                                 random_state=RANDOM_STATE)
        red.fit(X_sc, y_tr)
        return np.abs(red.coef_).mean(axis=0), scaler
    if rname == "EN":
        red = LogisticRegression(max_iter=2000, solver="saga", penalty="elasticnet",
                                 l1_ratio=0.5, random_state=RANDOM_STATE)
        red.fit(X_sc, y_tr)
        return np.abs(red.coef_).mean(axis=0), scaler
    raise ValueError(rname)


def make_classifiers():
    """name -> (template, needs_scale). Fresh templates (clone before fit)."""
    return {
        "RF": (RandomForestClassifier(n_estimators=100, max_depth=10,
                                      random_state=RANDOM_STATE, n_jobs=-1), False),
        "XGB": (XGBClassifier(n_estimators=200, max_depth=6, learning_rate=0.1,
                              subsample=0.8, random_state=RANDOM_STATE,
                              n_jobs=-1, verbosity=0), False),
        "LR": (LogisticRegression(max_iter=1000, random_state=RANDOM_STATE), True),
        "SVM_RBF": (SVC(kernel="rbf", probability=True, random_state=RANDOM_STATE), True),
        "SVM_Lin": (SVC(kernel="linear", probability=True, random_state=RANDOM_STATE), True),
        "GBM": (GradientBoostingClassifier(n_estimators=100, max_depth=5,
                                           random_state=RANDOM_STATE), False),
    }


# --------------------------------------------------------------------------
# Data + pool building
# --------------------------------------------------------------------------

class Data:
    def __init__(self):
        cfg, PATHS, POP, HARD, SIT, ML = setup()
        self.PATHS, self.SIT = PATHS, SIT
        d = np.load(str(PATHS.GENOTYPE_MATRIX))
        self.G = d["G"]
        self.samples = np.array(d["samples"].astype(str))
        with open(str(PATHS.GENOTYPE_MATRIX_COLS)) as fh:
            self.cols = [l.strip() for l in fh]
        self.snp_index = {s: i for i, s in enumerate(self.cols)}
        sm = pd.read_csv(str(PATHS.SAMPLES_CSV), header=None,
                         names=["sample", "pop", "super_pop"])
        lut = dict(zip(sm["sample"], sm["pop"]))
        self.y_str = np.array([lut[s] for s in self.samples])
        self.le = LabelEncoder().fit(self.y_str)
        self.y = self.le.transform(self.y_str)
        self.pfile = str(PATHS.PLINK_LD_PRUNED)
        self.fst_top_n = SIT.FST_TOP_N

    def build_pools(self, tr_idx, which=PANELS):
        """Return {panel: column-index array} built from training rows tr_idx only."""
        G, y_str = self.G, self.y_str
        pools = {}
        if "stat" in which:
            pools["stat"] = stat_pool_indices(G[tr_idx], y_str[tr_idx], top_n=500)
        fst_idx = None
        if "FST" in which or "fst_stat" in which:
            fst_idx = fst_pool_for_samples(
                self.samples[tr_idx], pfile=self.pfile, snp_col_index=self.snp_index,
                psam_ids=self.samples, top_n=self.fst_top_n)
        if "FST" in which:
            pools["FST"] = fst_idx
        if "fst_stat" in which:
            pools["fst_stat"] = fst_stat_pool_indices(G[tr_idx], y_str[tr_idx], fst_idx, self.cols)
        return pools


def _score(clf_tmpl, needs_scale, X_tr, y_tr, X_te, y_te):
    clf = clone(clf_tmpl)
    clf.fit(X_tr, y_tr)
    y_pred = clf.predict(X_te)
    rauc = None
    if hasattr(clf, "predict_proba"):
        try:
            rauc = roc_auc_score(y_te, clf.predict_proba(X_te),
                                 multi_class="ovr", average="macro")
        except Exception:
            pass
    return (accuracy_score(y_te, y_pred),
            f1_score(y_te, y_pred, average="weighted"),
            matthews_corrcoef(y_te, y_pred), rauc)


def _reduce(reductor, X_tr_pool, y_tr, N):
    """Fit reductor on the pool matrix; return top-N local indices + scaler."""
    imp, scaler = fit_reductor(reductor, X_tr_pool, y_tr)
    top = np.argsort(imp)[::-1][:N]
    return top, scaler


# --------------------------------------------------------------------------
# Tier B — isolate the leak (baseline config, in-fold pools)
# --------------------------------------------------------------------------

def run_tier_b(data: Data, baseline_dir: Path, ns=PANEL_SIZES):
    s1 = pd.read_csv(baseline_dir / "stage1_best.csv").set_index("n_snps")
    s2 = pd.read_csv(baseline_dir / "stage2_best.csv").set_index("n_snps")
    clsf = make_classifiers()
    skf = StratifiedKFold(n_splits=N_CV, shuffle=True, random_state=RANDOM_STATE)
    G, y = data.G, data.y

    # cache reductor importance per (fold, panel, reductor); pools per (fold, panel)
    folds = list(skf.split(G, y))
    pool_cache, red_cache = {}, {}
    needed_panels = sorted({s1.loc[N, "panel"] for N in ns})
    print(f"Tier B: panels needed = {needed_panels}")

    for fi, (tr, te) in enumerate(folds):
        t0 = time.time()
        pool_cache[fi] = data.build_pools(tr, which=needed_panels)
        print(f"  fold {fi+1}/{N_CV}: pools built in {time.time()-t0:.1f}s "
              f"({ {k: len(v) for k, v in pool_cache[fi].items()} })")

    records = []
    for N in ns:
        panel = s1.loc[N, "panel"]
        reductor = s1.loc[N, "reductor"]
        classifier = s2.loc[N, "classifier"]
        clf_tmpl, needs_scale = clsf[classifier]
        accs, f1s, mccs, raucs = [], [], [], []
        for fi, (tr, te) in enumerate(folds):
            pcols = pool_cache[fi][panel]
            X_tr_pool = G[tr][:, pcols].astype(np.float32)
            X_te_pool = G[te][:, pcols].astype(np.float32)
            key = (fi, panel, reductor)
            if key not in red_cache:
                red_cache[key] = fit_reductor(reductor, X_tr_pool, data.y[tr])
            imp, scaler = red_cache[key]
            top = np.argsort(imp)[::-1][:N]
            if needs_scale:
                X_tr = scaler.transform(X_tr_pool)[:, top]
                X_te = scaler.transform(X_te_pool)[:, top]
            else:
                X_tr = X_tr_pool[:, top]
                X_te = X_te_pool[:, top]
            a, f, m, r = _score(clf_tmpl, needs_scale, X_tr, data.y[tr], X_te, data.y[te])
            accs.append(a); f1s.append(f); mccs.append(m)
            if r is not None:
                raucs.append(r)
        records.append({
            "n_snps": N, "panel": panel, "reductor": reductor, "classifier": classifier,
            "acc_mean": np.mean(accs), "acc_std": np.std(accs),
            "f1_mean": np.mean(f1s), "mcc_mean": np.mean(mccs),
            "roc_auc_mean": np.mean(raucs) if raucs else np.nan,
        })
        print(f"  N={N:3d}  {panel}+{reductor}+{classifier}  "
              f"nested acc={np.mean(accs):.4f} ± {np.std(accs):.4f}")
    return pd.DataFrame(records)


# --------------------------------------------------------------------------
# Tier A — fully nested (inner-CV config selection)
# --------------------------------------------------------------------------

def run_tier_a(data: Data, ns, classifiers_subset=None):
    """Nested CV: pools on outer-train; config picked by inner CV; scored on outer-test.

    Inner CV reuses the outer-train pools (built once per outer fold) — this keeps
    the outer-test fold fully independent of pool/reductor/config/classifier, and
    only introduces negligible optimism into *config selection* (documented).
    """
    clsf = make_classifiers()
    if classifiers_subset:
        clsf = {k: v for k, v in clsf.items() if k in classifiers_subset}
    skf = StratifiedKFold(n_splits=N_CV, shuffle=True, random_state=RANDOM_STATE)
    inner_skf = StratifiedKFold(n_splits=N_CV, shuffle=True, random_state=RANDOM_STATE)
    G, y = data.G, data.y
    rf_eval = RandomForestClassifier(n_estimators=100, random_state=RANDOM_STATE, n_jobs=-1)

    records = []
    for fo, (tr, te) in enumerate(skf.split(G, y)):
        t0 = time.time()
        pools = data.build_pools(tr, which=PANELS)
        y_tr, y_te = y[tr], y[te]
        # Pre-extract pool matrices for outer-train / outer-test
        Xtr = {p: G[tr][:, c].astype(np.float32) for p, c in pools.items()}
        Xte = {p: G[te][:, c].astype(np.float32) for p, c in pools.items()}
        inner = list(inner_skf.split(Xtr["stat"], y_tr))
        print(f"[outer {fo+1}/{N_CV}] pools {{"
              f"{', '.join(f'{k}:{len(v)}' for k,v in pools.items())}}} "
              f"in {time.time()-t0:.1f}s")

        # cache reductor importances per (inner fold, panel, reductor) — N-independent
        red_imp = {}
        for ii, (itr, ival) in enumerate(inner):
            for panel in PANELS:
                Xp = Xtr[panel][itr]
                for rname in REDUCTORS:
                    red_imp[(ii, panel, rname)] = fit_reductor(rname, Xp, y_tr[itr])

        for N in ns:
            # ---- Stage-1-like: best (panel, reductor) by inner RF-eval acc ----
            best_pr, best_acc = None, -1
            for panel in PANELS:
                for rname in REDUCTORS:
                    accs = []
                    for ii, (itr, ival) in enumerate(inner):
                        imp, scaler = red_imp[(ii, panel, rname)]
                        top = np.argsort(imp)[::-1][:N]
                        Xa = Xtr[panel][itr][:, top]
                        Xb = Xtr[panel][ival][:, top]
                        e = clone(rf_eval); e.fit(Xa, y_tr[itr])
                        accs.append(e.score(Xb, y_tr[ival]))
                    if np.mean(accs) > best_acc:
                        best_acc, best_pr = np.mean(accs), (panel, rname)
            panel, rname = best_pr

            # ---- Stage-2-like: best classifier for chosen (panel, reductor) ----
            best_clf, best_cacc = None, -1
            for cname, (ctmpl, needs_scale) in clsf.items():
                accs = []
                for ii, (itr, ival) in enumerate(inner):
                    imp, scaler = red_imp[(ii, panel, rname)]
                    top = np.argsort(imp)[::-1][:N]
                    if needs_scale:
                        Xa = scaler.transform(Xtr[panel][itr])[:, top]
                        Xb = scaler.transform(Xtr[panel][ival])[:, top]
                    else:
                        Xa = Xtr[panel][itr][:, top]
                        Xb = Xtr[panel][ival][:, top]
                    a, *_ = _score(ctmpl, needs_scale, Xa, y_tr[itr], Xb, y_tr[ival])
                    accs.append(a)
                if np.mean(accs) > best_cacc:
                    best_cacc, best_clf = np.mean(accs), cname
            cname = best_clf
            ctmpl, needs_scale = clsf[cname]

            # ---- Refit on full outer-train, score untouched outer-test ----
            imp, scaler = fit_reductor(rname, Xtr[panel], y_tr)
            top = np.argsort(imp)[::-1][:N]
            if needs_scale:
                Xa = scaler.transform(Xtr[panel])[:, top]
                Xb = scaler.transform(Xte[panel])[:, top]
            else:
                Xa = Xtr[panel][:, top]
                Xb = Xte[panel][:, top]
            a, f, m, r = _score(ctmpl, needs_scale, Xa, y_tr, Xb, y_te)
            records.append({
                "outer_fold": fo, "n_snps": N, "panel": panel, "reductor": rname,
                "classifier": cname, "acc": a, "f1": f, "mcc": m, "roc_auc": r,
            })
            print(f"    N={N:3d}  {panel}+{rname}+{cname}  outer-test acc={a:.4f}")
    return pd.DataFrame(records)


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", choices=["B", "A", "both"], default="B")
    ap.add_argument("--ns", default="", help="comma list of N (default: full curve / committed for A)")
    args = ap.parse_args()

    data = Data()
    out_dir = Path(str(data.PATHS.outputs_dir("self_evaluation/08b_nested_cv_sweep")))
    out_dir.mkdir(parents=True, exist_ok=True)
    baseline_dir = (Path(__file__).resolve().parent.parent
                    / "results_archive/baseline_leaky_v1/self_evaluation/08_unified_panel_sweep")

    if args.tier in ("B", "both"):
        ns = [int(x) for x in args.ns.split(",")] if args.ns else PANEL_SIZES
        print("=" * 70, "\nTIER B — isolate leak (baseline config, in-fold pools)\n", "=" * 70)
        df = run_tier_b(data, baseline_dir, ns)
        p = out_dir / "nested_tierB_results.csv"
        df.to_csv(p, index=False)
        print(f"Saved: {p}")

    if args.tier in ("A", "both"):
        ns = [int(x) for x in args.ns.split(",")] if args.ns else [35, 50, 70]
        print("=" * 70, "\nTIER A — fully nested (inner-CV config selection)\n", "=" * 70)
        df = run_tier_a(data, ns)
        agg = (df.groupby(["n_snps"])
                 .agg(acc_mean=("acc", "mean"), acc_std=("acc", "std"),
                      f1_mean=("f1", "mean"), mcc_mean=("mcc", "mean"),
                      roc_auc_mean=("roc_auc", "mean"))
                 .reset_index())
        df.to_csv(out_dir / "nested_tierA_folds.csv", index=False)
        agg.to_csv(out_dir / "nested_tierA_results.csv", index=False)
        print(f"Saved: {out_dir / 'nested_tierA_results.csv'}")
        print(agg.to_string(index=False))


if __name__ == "__main__":
    main()
