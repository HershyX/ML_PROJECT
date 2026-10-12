"""Precompute PCA + k-means for the TCGA RNA-seq matrix and export data.js for the frontend."""
import json
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

df = pd.read_csv("../RNA_seq_1.csv", encoding="utf-8-sig", index_col=0)
genes = df.index.tolist()
samples = df.columns.tolist()
X = df.T.to_numpy(dtype=float)  # samples x genes
print("samples x genes", X.shape, "missing", int(np.isnan(X).sum()))

n_missing = int(np.isnan(X).sum())
mu, sd = np.nanmean(X, 0), np.nanstd(X, axis=0)
sd[sd == 0] = 1
Z = np.nan_to_num((X - mu) / sd)  # missing -> gene mean
X = np.where(np.isnan(X), mu, X)

pca = PCA(n_components=10, random_state=0)
P = pca.fit_transform(Z)

out = {
    "nMissing": n_missing,
    "samples": samples,
    "genes": genes,
    "geneMean": np.round(mu, 3).tolist(),
    "geneVar": np.round(X.var(0), 4).tolist(),
    "pc": np.round(P[:, :5], 3).tolist(),
    "explained": np.round(pca.explained_variance_ratio_ * 100, 2).tolist(),
    "sampleSd": np.round(X.std(1), 4).tolist(),
    "k": {},
}

order = np.argsort(-X.var(0))
top_idx = order[:40]
out["topVarGenes"] = [{"gene": genes[i], "var": round(float(X.var(0)[i]), 3)} for i in order[:20]]
out["heatGenes"] = [genes[i] for i in top_idx]
out["heat"] = np.round(Z[:, top_idx], 2).T.tolist()  # gene x sample, z-scored

for k in range(2, 7):
    km = KMeans(n_clusters=k, n_init=20, random_state=0).fit(P)
    lab = km.labels_
    # relabel clusters by size, largest first
    rank = {c: r for r, c in enumerate(np.argsort(-np.bincount(lab)))}
    lab = np.array([rank[c] for c in lab])
    cmeans = np.array([Z[lab == c].mean(0) for c in range(k)])  # k x genes (z-scored)
    markers = []
    for c in range(k):
        diff = cmeans[c] - np.delete(cmeans, c, 0).mean(0)
        top = np.argsort(-diff)[:8]
        markers.append([{"gene": genes[i], "z": round(float(cmeans[c][i]), 2)} for i in top])
    out["k"][str(k)] = {
        "labels": lab.tolist(),
        "sizes": np.bincount(lab).tolist(),
        "silhouette": round(float(silhouette_score(P, lab)), 4),
        "inertia": round(float(km.inertia_), 1),
        "markers": markers,
        "geneClusterZ": np.round(cmeans, 2).T.tolist(),  # genes x k
    }
    print(k, out["k"][str(k)]["sizes"], out["k"][str(k)]["silhouette"])

with open("data.js", "w") as f:
    f.write("window.TCGA = " + json.dumps(out, separators=(",", ":")) + ";")
print("wrote data.js")
