"""Fetch rsIDs for the committed-panel marker table via Ensembl's GRCh37
VEP REST API (single batch POST, <=200 variants per call — well within
our 94 unique SNPs). Coordinates in committed_panel_markers.csv are b37,
so the GRCh37-specific endpoint (grch37.rest.ensembl.org) is used, not the
default GRCh38 one -- using the wrong build would silently mismatch.

Writes: reports/figures/committed_panel_markers.csv (adds an rsID column,
overwriting in place)
"""

import json
import time
import urllib.request
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
CSV = REPO / "reports/figures/committed_panel_markers.csv"
ENDPOINT = "https://grch37.rest.ensembl.org/vep/human/region"

df = pd.read_csv(CSV)
uniq = df[["chr", "pos", "ref", "alt"]].drop_duplicates().reset_index(drop=True)
print(f"{len(uniq)} unique SNPs to look up")

variants = [f"{r.chr} {r.pos} . {r.ref} {r.alt} . . ." for r in uniq.itertuples()]

req = urllib.request.Request(
    ENDPOINT,
    data=json.dumps({"variants": variants}).encode(),
    headers={"Content-Type": "application/json", "Accept": "application/json"},
    method="POST",
)
with urllib.request.urlopen(req, timeout=60) as resp:
    results = json.loads(resp.read())

assert len(results) == len(uniq), f"expected {len(uniq)} results, got {len(results)}"

rsid_map = {}
n_found = 0
for row, res in zip(uniq.itertuples(), results):
    key = (row.chr, row.pos)
    rsid = None
    for cv in res.get("colocated_variants", []):
        if cv.get("id", "").startswith("rs"):
            rsid = cv["id"]
            break
    rsid_map[key] = rsid
    if rsid:
        n_found += 1

print(f"rsIDs found: {n_found} / {len(uniq)}")

df["rsID"] = [rsid_map.get((r.chr, r.pos)) for r in df.itertuples()]
df.to_csv(CSV, index=False)
print(f"Updated {CSV}")
print(df[["N", "rank", "chr", "pos", "rsID"]].head(10).to_string(index=False))
