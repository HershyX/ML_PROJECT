"""Draw the four poster graphs from real pipeline results. Run: python posters/make_charts.py"""
import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "ui"))
import pipeline as pl  # noqa: E402

OUT = HERE / "assets"
BLUE, ORANGE, GREEN, RED, PURPLE = "#1f6fe0", "#f59e0b", "#22a55a", "#e5383b", "#7c3aed"
plt.rcParams.update({
    "font.family": "Segoe UI", "font.size": 9, "axes.titlesize": 10.5, "axes.titleweight": "bold",
    "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#9aa3b2",
    "axes.grid": True, "grid.color": "#e6e9ef", "grid.linewidth": 0.8, "axes.axisbelow": True,
    "xtick.color": "#4b5563", "ytick.color": "#4b5563", "axes.labelcolor": "#374151",
    "legend.frameon": False, "svg.fonttype": "none",
})


def save(fig, name):
    fig.tight_layout(pad=0.6)
    fig.savefig(OUT / f"{name}.svg")
    plt.close(fig)


def silhouette_by_k():
    configs = json.loads((ROOT / "galaxy" / "data.json").read_text(encoding="utf-8"))["configs"]
    ks = np.arange(2, 7)
    km = [configs[f"kmeans-{k}"]["metrics"]["silhouette"] for k in ks]
    hc = [configs[f"hierarchical-{k}"]["metrics"]["silhouette"] for k in ks]
    fig, ax = plt.subplots(figsize=(3.6, 2.5))
    ax.bar(ks - 0.19, km, 0.38, color=BLUE, label="K-Means")
    ax.bar(ks + 0.19, hc, 0.38, color=ORANGE, label="Hierarchical")
    ax.set(title="Silhouette Score by Number of Groups", xlabel="Number of groups (k)", ylabel="Silhouette")
    ax.set_xticks(ks)
    ax.legend(fontsize=8)
    save(fig, "graph_silhouette")


def explained_variance(result):
    evr = result["explained_variance_ratio"]
    comps = np.arange(1, len(evr) + 1)
    fig, ax = plt.subplots(figsize=(3.6, 2.5))
    ax.bar(comps, evr * 100, color=BLUE, alpha=0.55, label="Each component")
    ax.plot(comps, np.cumsum(evr) * 100, color=GREEN, lw=2, label="Cumulative")
    ax.set(title="PCA Explained Variance", xlabel="Principal component", ylabel="Variance explained (%)", ylim=(0, 80))
    ax.legend(fontsize=8, loc="center right")
    save(fig, "graph_variance")


def unseen_confidence(result):
    clf = pl.PatientClassifier(result)
    X = pl.load_expression(ROOT / "test_data" / "test_set.csv")
    labels = [clf.classify(X.loc[p])["confidence_label"] for p in X.index]
    names = ["Clear match", "Moderate match", "Borderline"]
    counts = [labels.count(n) for n in names]
    fig, ax = plt.subplots(figsize=(3.6, 2.5))
    bars = ax.bar(names, counts, color=[GREEN, ORANGE, RED], width=0.6)
    ax.bar_label(bars, padding=2, fontsize=9, fontweight="bold")
    ax.set(title=f"Unseen Test Patients (n = {len(labels)})", ylabel="Patients", ylim=(0, max(counts) * 1.2))
    save(fig, "graph_unseen")


def execution_time(X):
    prep = pl.Preprocessor().fit(X)
    Z = prep.transform(X)
    order = np.argsort(-X[prep.genes].var().to_numpy())
    sizes = [500, 1000, 2000, 3000, len(order)]

    def best_of(fn, runs=3):
        times = []
        for _ in range(runs):
            t = time.perf_counter()
            fn()
            times.append(time.perf_counter() - t)
        return min(times)

    def cluster_and_score(space):
        labels = pl.cluster_labels(space, "kmeans", 4)
        silhouette_score(space, labels)

    raw, reduced = [], []
    for n in sizes:
        sub = Z[:, order[:n]]
        raw.append(best_of(lambda: cluster_and_score(sub)))
        reduced.append(best_of(lambda: cluster_and_score(PCA(n_components=50, random_state=pl.SEED).fit_transform(sub))))
    fig, ax = plt.subplots(figsize=(3.6, 2.5))
    ax.plot(sizes, raw, marker="o", color=RED, lw=2, label="All genes")
    ax.plot(sizes, reduced, marker="o", color=PURPLE, lw=2, label="PCA first")
    ax.set(title="Execution Time vs Number of Genes", xlabel="Genes used", ylabel="Time (seconds)")
    ax.legend(fontsize=8)
    save(fig, "graph_time")
    print("time", [round(v, 3) for v in raw], [round(v, 3) for v in reduced])


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    X = pl.load_expression(ROOT / "RNA_seq_1.csv")
    result = pl.run_pipeline(X)
    silhouette_by_k()
    explained_variance(result)
    unseen_confidence(result)
    execution_time(X.loc[result["patient_ids"]])
    print("charts written to", OUT)
