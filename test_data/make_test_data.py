"""Cut small testing files out of the real TCGA data (RNA_seq_1.csv). No values are made up.

Writes, in the same layout as the original (genes as rows, patients as columns):
- test_set.csv: the patients both apps hold back (same seed and 20% split), all genes.
- small_sample.csv: 60 patients and the 300 most variable genes, for a quick upload test.
- extra_patient.csv: one of the held-back patients on its own, to try importing a patient.

Run with: python test_data/make_test_data.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SEED = 42
HOLDOUT = 0.2

raw = pd.read_csv(HERE.parent / "RNA_seq_1.csv", encoding="utf-8-sig")
gene_col = raw.columns[0]
patients = raw.columns[1:]

held = np.random.default_rng(SEED).random(len(patients)) < HOLDOUT  # matches ui/pipeline.py
raw[[gene_col, *patients[held]]].to_csv(HERE / "test_set.csv", index=False)

rng = np.random.default_rng(SEED + 1)
picked = sorted(rng.choice(len(patients), size=60, replace=False))
top_genes = raw[patients].var(axis=1).nlargest(300).index.sort_values()
raw.loc[top_genes, [gene_col, *patients[picked]]].to_csv(HERE / "small_sample.csv", index=False)

# One unseen patient on its own, to try the import in both apps
raw[[gene_col, patients[held][0]]].to_csv(HERE / "extra_patient.csv", index=False)

for name in ("test_set.csv", "small_sample.csv", "extra_patient.csv"):
    df = pd.read_csv(HERE / name)
    print(f"{name}: {df.shape[0]} genes, {df.shape[1] - 1} patients, {(HERE / name).stat().st_size / 1e6:.2f} MB")
