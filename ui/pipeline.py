"""Model code for the genomic subtype clustering app.

Pipeline from the project abstract: standardize, PCA, t-SNE, then K-Means or
hierarchical clustering, scored on cluster quality and run time.
"""
from __future__ import annotations

import string
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.metrics import (adjusted_rand_score, calinski_harabasz_score,
                             davies_bouldin_score, silhouette_score)
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler

SEED = 42
GROUP_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]


def group_name(i: int) -> str:
    return f"Group {string.ascii_uppercase[i]}"


def load_expression(source) -> pd.DataFrame:
    """Read a genes x patients CSV and return a patients x genes table."""
    df = pd.read_csv(source, encoding="utf-8-sig")
    if df.shape[1] < 2:
        raise ValueError("Expected a gene name column followed by at least one patient column.")
    gene_col = df.columns[0]
    values = df.drop(columns=gene_col).apply(pd.to_numeric, errors="coerce")
    if values.isna().mean().mean() > 0.5:
        raise ValueError("Most values are not numbers. Check that genes are rows and patients are columns.")
    values.index = df[gene_col].astype(str).str.strip()
    values = values.groupby(level=0, sort=False).mean()  # average repeated gene names
    X = values.T
    X.index = X.index.astype(str)
    X.index.name = "patient_id"
    X.columns.name = None
    return X


@dataclass
class Preprocessor:
    """Fill blanks with each gene's median, drop flat genes, standardize."""
    genes: list[str] = field(default_factory=list)
    medians: pd.Series | None = None
    scaler: StandardScaler | None = None

    def fit(self, X: pd.DataFrame) -> "Preprocessor":
        medians = X.median()
        filled = _fill(X.to_numpy(dtype=float), medians.to_numpy())
        keep = np.nan_to_num(filled.var(axis=0)) > 1e-8  # all-blank genes come out as NaN
        self.genes = list(X.columns[keep])
        self.medians = medians[self.genes]
        self.scaler = StandardScaler().fit(filled[:, keep])
        return self

    def transform(self, X: pd.DataFrame) -> np.ndarray:
        aligned = X.reindex(columns=self.genes).to_numpy(dtype=float)
        return self.scaler.transform(_fill(aligned, self.medians.to_numpy()))


def _fill(values: np.ndarray, medians: np.ndarray) -> np.ndarray:
    # numpy instead of DataFrame.fillna, which goes column by column and is slow on 4,571 genes
    return np.where(np.isnan(values), medians, values)


def cluster_labels(space: np.ndarray, method: str, k: int) -> np.ndarray:
    if method == "kmeans":
        model = KMeans(n_clusters=k, n_init=10, random_state=SEED)
    else:
        model = AgglomerativeClustering(n_clusters=k, linkage="ward")
    return model.fit_predict(space)


def cluster_scores(space: np.ndarray, labels: np.ndarray) -> dict:
    return {
        "silhouette": float(silhouette_score(space, labels)),
        "davies_bouldin": float(davies_bouldin_score(space, labels)),
        "calinski_harabasz": float(calinski_harabasz_score(space, labels)),
    }


def order_by_size(labels: np.ndarray) -> np.ndarray:
    """Renumber groups so 0 is the largest."""
    order = np.argsort(-np.bincount(labels), kind="stable")
    remap = np.empty_like(order)
    remap[order] = np.arange(len(order))
    return remap[labels]


def run_pipeline(X: pd.DataFrame, n_pca: int = 50, perplexity: float = 30,
                 method: str = "kmeans", k: int = 4, cluster_on: str = "pca",
                 holdout: float = 0.2) -> dict:
    """Fit the whole pipeline on all patients except a random `holdout` share,
    which is kept back as unseen patients for the classifier."""
    rng = np.random.default_rng(SEED)
    held = rng.random(len(X)) < holdout
    holdout_ids = list(X.index[held])
    X = X[~held]
    n = len(X)
    if n <= k:
        raise ValueError(f"Need more than {k} patients to make {k} groups.")
    timings = {}

    t = time.perf_counter()
    prep = Preprocessor().fit(X)
    Z = prep.transform(X)
    timings["preprocess"] = time.perf_counter() - t

    n_pca = max(2, min(n_pca, min(Z.shape) - 1))
    t = time.perf_counter()
    pca = PCA(n_components=n_pca, random_state=SEED)
    pca_coords = pca.fit_transform(Z)
    timings["pca"] = time.perf_counter() - t

    t = time.perf_counter()
    tsne = TSNE(n_components=2, perplexity=min(perplexity, n - 1), init="pca", random_state=SEED)
    tsne_coords = tsne.fit_transform(pca_coords)
    timings["tsne"] = time.perf_counter() - t

    space = {"pca": pca_coords, "tsne": tsne_coords, "raw": Z}[cluster_on]
    t = time.perf_counter()
    labels = order_by_size(cluster_labels(space, method, k))
    timings["clustering"] = time.perf_counter() - t
    timings["total"] = sum(timings.values())

    return {
        "preprocessor": prep,
        "pca": pca,
        "Z": Z,
        "space": space,
        "pca_coords": pca_coords,
        "tsne_coords": tsne_coords,
        "explained_variance_ratio": pca.explained_variance_ratio_,
        "labels": labels,
        "metrics": cluster_scores(space, labels),
        "timings": timings,
        "patient_ids": list(X.index),
        "holdout_ids": holdout_ids,
        "genes": prep.genes,
        "settings": dict(n_pca=n_pca, perplexity=perplexity, method=method, k=k,
                         cluster_on=cluster_on, holdout=holdout),
    }


def sweep_k(space: np.ndarray, method: str, ks=range(2, 9)) -> pd.DataFrame:
    """Cluster quality for each number of groups."""
    rows = [{"k": k, **cluster_scores(space, cluster_labels(space, method, k))} for k in ks if k < len(space)]
    return pd.DataFrame(rows)


def speed_comparison(Z: np.ndarray, n_pca: int, k: int, method: str,
                     perplexity: float = 30) -> tuple[pd.DataFrame, float]:
    """Time each step on all genes against PCA first.

    The PCA columns include the time to fit PCA itself. Also returns the
    adjusted Rand index between the two groupings, i.e. whether reducing
    first changes who ends up together.
    """
    def timed(fn):
        t = time.perf_counter()
        out = fn()
        return out, time.perf_counter() - t

    reduced, pca_time = timed(lambda: PCA(n_components=n_pca, random_state=SEED).fit_transform(Z))

    def cluster_and_score(space):
        labels = cluster_labels(space, method, k)
        return labels, silhouette_score(space, labels)

    (raw_labels, raw_sil), raw_cluster = timed(lambda: cluster_and_score(Z))
    (pca_labels, pca_sil), pca_cluster = timed(lambda: cluster_and_score(reduced))

    perplexity = min(perplexity, len(Z) - 1)
    _, raw_tsne = timed(lambda: TSNE(perplexity=perplexity, init="pca", random_state=SEED).fit_transform(Z))
    _, pca_tsne = timed(lambda: TSNE(perplexity=perplexity, init="pca", random_state=SEED).fit_transform(reduced))

    table = pd.DataFrame([
        {"step": "Clustering and scoring", "all_genes_seconds": raw_cluster,
         "pca_first_seconds": pca_time + pca_cluster, "all_genes_silhouette": raw_sil,
         "pca_first_silhouette": pca_sil},
        {"step": "t-SNE map", "all_genes_seconds": raw_tsne,
         "pca_first_seconds": pca_time + pca_tsne},
    ])
    table["speedup"] = table["all_genes_seconds"] / table["pca_first_seconds"]
    return table, float(adjusted_rand_score(raw_labels, pca_labels))


def marker_genes(Z: np.ndarray, genes: list[str], labels: np.ndarray, top: int = 10) -> dict[int, pd.DataFrame]:
    """Genes most raised in each group compared with everyone else, in standard deviations."""
    genes = np.asarray(genes)
    out = {}
    for g in np.unique(labels):
        inside = labels == g
        diff = Z[inside].mean(axis=0) - Z[~inside].mean(axis=0)
        idx = np.argsort(-diff)[:top]
        out[int(g)] = pd.DataFrame({"gene": genes[idx], "difference": diff[idx].round(2)})
    return out


class PatientClassifier:
    """Place a new patient into one of the groups found by the pipeline.

    t-SNE cannot place new points, so the patient is projected with the fitted
    preprocessing and PCA, then the 15 nearest patients vote, weighted by
    closeness. The map position is the average t-SNE position of the 10
    nearest patients.
    """

    def __init__(self, result: dict, n_votes: int = 15, n_map: int = 10):
        self.r = result
        self.ids = np.asarray(result["patient_ids"])
        self.labels = result["labels"]
        self.k = int(self.labels.max()) + 1
        self.n_votes = min(n_votes, len(self.ids))
        self.n_map = min(n_map, self.n_votes)
        self.nn = NearestNeighbors().fit(result["pca_coords"])
        self.markers = marker_genes(result["Z"], result["genes"], self.labels, top=20)
        self.gene_index = {g: i for i, g in enumerate(result["genes"])}

    def classify(self, profile: pd.Series) -> dict:
        prep = self.r["preprocessor"]
        profile = pd.to_numeric(profile, errors="coerce")
        profile.index = profile.index.astype(str).str.strip()
        known = profile.dropna().groupby(level=0).mean()
        overlap = known.index.intersection(prep.genes)
        if len(overlap) < 10:
            raise ValueError(f"Only {len(overlap)} of the uploaded genes match the dataset. Check the gene names.")

        z = prep.transform(known.to_frame().T)
        point = self.r["pca"].transform(z)

        dist, idx = self.nn.kneighbors(point, n_neighbors=self.n_votes)
        dist, idx = dist[0], idx[0]
        weights = 1.0 / (dist + 1e-9)
        proba = np.bincount(self.labels[idx], weights=weights, minlength=self.k) / weights.sum()
        group = int(np.argmax(proba))
        confidence = float(proba[group])
        if confidence >= 0.8:
            confidence_label = "Clear match"
        elif confidence >= 0.5:
            confidence_label = "Moderate match"
        else:
            confidence_label = "Borderline"

        map_idx = idx[:self.n_map]
        neighbors = pd.DataFrame({
            "patient_id": self.ids[idx[:5]],
            "group": [group_name(g) for g in self.labels[idx[:5]]],
            "distance": dist[:5].round(2),
        })

        markers = self.markers[group]["gene"]
        patient_z = np.array([z[0, self.gene_index[g]] for g in markers])
        direction = np.where(patient_z > 0.5, "Higher than average",
                             np.where(patient_z < -0.5, "Lower than average", "About average"))
        driving = (pd.DataFrame({"gene": markers.to_numpy(), "patient_z": patient_z.round(2), "direction": direction})
                   .sort_values("patient_z", ascending=False).head(8).reset_index(drop=True))

        return {
            "group": group,
            "confidence": confidence,
            "confidence_label": confidence_label,
            "probabilities": proba,
            "neighbors": neighbors,
            "tsne_xy": self.r["tsne_coords"][map_idx].mean(axis=0),
            "pca_xy": point[0, :2],
            "driving_genes": driving,
            "n_missing_genes": len(prep.genes) - len(overlap),
            "n_extra_genes": len(known.index.difference(prep.genes)),
        }
