"""Precompute everything the 3D patient galaxy shows and write data.json.

Uses the same model code as the Streamlit app (ui/pipeline.py). Run with:
python galaxy/build_data.py
"""
import base64
import json
import sys
from pathlib import Path

import numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.neighbors import NearestNeighbors

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "ui"))
import pipeline as pl  # noqa: E402

N_PCA = 50
HOLDOUT = 0.2
KS = range(2, 7)
METHODS = ("kmeans", "hierarchical")
N_VOTES = 15
N_MAP = 10
N_VARIABLE_GENES = 400  # kept for colour by gene, on top of every marker gene
Z_SCALE = 16  # standardized values are stored as int8 of z * Z_SCALE


def normalise(coords, center=None, scale=None):
    """Center the cloud and scale it so 95% of patients sit within radius 60."""
    if center is None:
        center = coords.mean(axis=0)
        scale = np.percentile(np.linalg.norm(coords - center, axis=1), 95) / 60
    return (coords - center) / scale, center, scale


def flat(a):
    return [round(float(v), 2) for v in np.ravel(a)]


def f32(a):
    """float32 array as base64, decoded in the browser into a Float32Array."""
    return base64.b64encode(np.ascontiguousarray(a, dtype=np.float32).tobytes()).decode()


def main():
    X = pl.load_expression(HERE.parent / "RNA_seq_1.csv")
    held = np.random.default_rng(pl.SEED).random(len(X)) < HOLDOUT  # same split as the Streamlit app
    train, test = X[~held], X[held]

    prep = pl.Preprocessor().fit(train)
    Z_train, Z_test = prep.transform(train), prep.transform(test)
    pca = PCA(n_components=N_PCA, random_state=pl.SEED).fit(Z_train)
    P_train, P_test = pca.transform(Z_train), pca.transform(Z_test)

    tsne = TSNE(n_components=3, perplexity=30, init="pca", random_state=pl.SEED).fit_transform(P_train)
    tsne, _, _ = normalise(tsne)
    pca3, center, scale = normalise(P_train[:, :3])
    pca3_test, _, _ = normalise(P_test[:, :3], center, scale)

    nn = NearestNeighbors().fit(P_train)
    own = nn.kneighbors(P_train, n_neighbors=6, return_distance=False)[:, 1:]  # drop the patient itself
    vote_dist, vote_idx = nn.kneighbors(P_test, n_neighbors=N_VOTES)
    tsne_test = tsne[vote_idx[:, :N_MAP]].mean(axis=1)  # t-SNE cannot place new points

    configs, marker_names = {}, set()
    for method in METHODS:
        for k in KS:
            labels = pl.order_by_size(pl.cluster_labels(P_train, method, k))
            markers = pl.marker_genes(Z_train, prep.genes, labels, top=8)
            configs[f"{method}-{k}"] = {
                "labels": labels.tolist(),
                "metrics": {key: round(v, 3) for key, v in pl.cluster_scores(P_train, labels).items()},
                "markers": [[[row.gene, float(row.difference)] for row in markers[g].itertuples()]
                            for g in range(k)],
            }
            marker_names.update(gene for table in markers.values() for gene in table["gene"])

    variance = train[prep.genes].var().sort_values(ascending=False)
    names = sorted(marker_names | set(variance.index[:N_VARIABLE_GENES]))
    cols = [prep.genes.index(g) for g in names]
    Z_all = np.vstack([Z_train, Z_test])[:, cols]  # training patients first, then unseen ones
    packed = np.clip(np.round(Z_all * Z_SCALE), -127, 127).astype(np.int8)

    data = {
        "meta": {
            "patients": len(X), "genes": X.shape[1], "train": len(train), "unseen": len(test),
            "n_pca": N_PCA, "pca3_variance": round(float(pca.explained_variance_ratio_[:3].sum()), 3),
        },
        "train": {"ids": list(train.index), "tsne": flat(tsne), "pca": flat(pca3), "neighbors": own.tolist()},
        "unseen": {"ids": list(test.index), "tsne": flat(tsne_test), "pca": flat(pca3_test),
                   "vote_idx": vote_idx.tolist(), "vote_dist": np.round(vote_dist, 2).tolist()},
        "configs": configs,
        "genes": {"names": names, "scale": Z_SCALE, "z": base64.b64encode(packed.tobytes()).decode()},
    }
    # Everything the page needs to place an imported patient the same way the Python model does
    model = {
        "genes": prep.genes, "k": N_PCA, "votes": N_VOTES, "map_neighbours": N_MAP,
        "median": f32(prep.medians.to_numpy()), "mean": f32(prep.scaler.mean_), "scale": f32(prep.scaler.scale_),
        "pca_mean": f32(pca.mean_), "components": f32(pca.components_), "train_pca": f32(P_train),
        "pca3_center": center.tolist(), "pca3_scale": float(scale),
    }
    (HERE / "model.json").write_text(json.dumps(model, separators=(",", ":")), encoding="utf-8")

    out = HERE / "data.json"
    out.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    print(f"wrote {out.name}: {out.stat().st_size / 1e6:.1f} MB, {len(train)} grouped, {len(test)} unseen, "
          f"{len(names)} genes, {len(configs)} settings")


if __name__ == "__main__":
    main()
