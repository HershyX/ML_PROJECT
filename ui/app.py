"""Streamlit app for the genomic subtype clustering project.

Run from the ML_PROJECT folder with: streamlit run ui/app.py
"""
import io
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

import pipeline as pl

DATA_PATH = Path(__file__).resolve().parent.parent / "RNA_seq_1.csv"
METHODS = {"K-Means": "kmeans", "Hierarchical (Ward)": "hierarchical"}
SPACES = {"PCA output (recommended)": "pca", "t-SNE map": "tsne", "All genes (no reduction)": "raw"}
FOOTER = "For research and teaching. Groups are statistical, not diagnoses."

st.set_page_config(page_title="Genomic Subtype Clustering", layout="wide")


@st.cache_data(show_spinner="Loading data...")
def load_builtin() -> pd.DataFrame:
    return pl.load_expression(DATA_PATH)


@st.cache_data(show_spinner="Loading data...")
def load_upload(data: bytes) -> pd.DataFrame:
    return pl.load_expression(io.BytesIO(data))


def scatter(ax, coords, labels, alpha=0.85, size=14, legend=True):
    for g in np.unique(labels):
        m = labels == g
        ax.scatter(coords[m, 0], coords[m, 1], s=size, alpha=alpha, linewidths=0,
                   color=pl.GROUP_COLORS[g], label=f"{pl.group_name(g)} ({m.sum()})")
    if legend:
        ax.legend(fontsize=8, frameon=False, markerscale=1.5)


def show(fig):
    fig.tight_layout()
    st.pyplot(fig)
    plt.close(fig)


# Sidebar

with st.sidebar:
    st.header("Settings")
    source = st.radio("Data", ["Built-in TCGA dataset", "Upload my own CSV"])
    X, data_key = None, source
    if source == "Built-in TCGA dataset":
        X = load_builtin()
    else:
        upload = st.file_uploader("Genes as rows, patients as columns", type="csv")
        if upload is not None:
            data_key = f"{upload.name}:{upload.size}"
            try:
                X = load_upload(upload.getvalue())
            except ValueError as e:
                st.error(str(e))
    if X is not None:
        st.caption(f"{X.shape[0]:,} patients, {X.shape[1]:,} genes")

    n_pca = st.slider("PCA components", 2, 100, 50)
    perplexity = st.slider("t-SNE perplexity", 5, 50, 30)
    method_label = st.selectbox("Grouping method", list(METHODS))
    k = st.slider("Number of groups", 2, 8, 4)
    space_label = st.selectbox("Cluster on", list(SPACES))
    holdout = st.slider("Patients held back for testing", 0, 40, 20, step=5, format="%d%%",
                        help="These patients are left out of the grouping and used as unseen patients "
                             "in the Classify a patient tab.")
    settings = dict(data=data_key, n_pca=n_pca, perplexity=perplexity, method=METHODS[method_label],
                    k=k, cluster_on=SPACES[space_label], holdout=holdout / 100)
    run = st.button("Run pipeline", type="primary", disabled=X is None)
    if "settings" in st.session_state and st.session_state.settings != settings and not run:
        st.caption("Settings changed. Run again to update.")

st.title("Genomic Subtype Clustering")
st.caption("Group patients by gene activity using PCA, t-SNE and clustering.")

if run:
    params = {key: value for key, value in settings.items() if key != "data"}
    with st.spinner("Running PCA, t-SNE and clustering..."):
        try:
            result = pl.run_pipeline(X, **params)
        except ValueError as e:
            st.error(str(e))
        else:
            st.session_state.update(result=result, settings=settings, X=X,
                                    classifier=pl.PatientClassifier(result))
            for key in ("sweep", "speed", "classified", "holdout_table"):
                st.session_state.pop(key, None)

result = st.session_state.get("result")
tab_run, tab_classify = st.tabs(["Run pipeline", "Classify a patient"])

# Tab 1: the pipeline

with tab_run:
    if result is None:
        st.info("Choose settings in the sidebar and press Run pipeline.")
    else:
        m, t, labels = result["metrics"], result["timings"], result["labels"]
        cols = st.columns(4)
        cols[0].metric("Silhouette", f"{m['silhouette']:.3f}",
                       help="From -1 to 1. Higher means better separated groups.")
        cols[1].metric("Davies-Bouldin", f"{m['davies_bouldin']:.2f}", help="Lower is better.")
        cols[2].metric("Calinski-Harabasz", f"{m['calinski_harabasz']:.1f}", help="Higher is better.")
        cols[3].metric("Total time", f"{t['total']:.1f} s",
                       help=f"Preprocess {t['preprocess']:.1f} s, PCA {t['pca']:.1f} s, "
                            f"t-SNE {t['tsne']:.1f} s, clustering {t['clustering']:.2f} s.")
        if m["silhouette"] >= 0.5:
            st.write("Groups are clearly separated.")
        elif m["silhouette"] >= 0.25:
            st.write("Groups are moderately separated.")
        else:
            st.write("Groups overlap a lot, so treat the boundaries as soft.")
        st.caption(f"Grouped {len(result['patient_ids'])} patients. "
                   f"{len(result['holdout_ids'])} were held back for testing.")

        evr = result["explained_variance_ratio"]
        left, right = st.columns(2)
        with left:
            fig, ax = plt.subplots(figsize=(5.5, 4.2))
            scatter(ax, result["pca_coords"], labels)
            ax.set_title("PCA map (PC1 vs PC2)", fontsize=11)
            ax.set_xlabel(f"PC1 ({evr[0]:.1%} of variation)")
            ax.set_ylabel(f"PC2 ({evr[1]:.1%} of variation)")
            show(fig)
        with right:
            fig, ax = plt.subplots(figsize=(5.5, 4.2))
            scatter(ax, result["tsne_coords"], labels)
            ax.set_title("t-SNE map", fontsize=11)
            ax.set_xticks([])
            ax.set_yticks([])
            show(fig)

        left, right = st.columns([3, 2])
        with left:
            st.subheader("Explained variance")
            cumulative = np.cumsum(evr)
            fig, ax = plt.subplots(figsize=(6, 3))
            comps = np.arange(1, len(evr) + 1)
            ax.bar(comps, evr * 100, color="#1f77b4", label="Each component")
            ax.plot(comps, cumulative * 100, color="#d62728", marker="o", markersize=3, label="Cumulative")
            ax.set_xlabel("Component")
            ax.set_ylabel("% of variation")
            ax.legend(fontsize=8, frameon=False)
            show(fig)
            note = f"The {len(evr)} components kept explain {cumulative[-1]:.0%} of the variation."
            for target in (0.8, 0.9):
                if cumulative[-1] >= target:
                    note += f" {target:.0%} is reached at {int(np.argmax(cumulative >= target)) + 1} components."
            st.caption(note)
        with right:
            st.subheader("Group sizes")
            counts = np.bincount(labels)
            st.dataframe(pd.DataFrame({
                "group": [pl.group_name(g) for g in range(len(counts))],
                "patients": counts,
                "percent": (counts / counts.sum() * 100).round(1),
            }), hide_index=True)

        st.subheader("Choosing the number of groups")
        if st.button("Compare 2 to 8 groups"):
            with st.spinner("Clustering for each number of groups..."):
                st.session_state.sweep = pl.sweep_k(result["space"], result["settings"]["method"])
        if "sweep" in st.session_state:
            sweep = st.session_state.sweep
            left, right = st.columns(2)
            for col, metric, title in ((left, "silhouette", "Silhouette (higher is better)"),
                                       (right, "davies_bouldin", "Davies-Bouldin (lower is better)")):
                with col:
                    fig, ax = plt.subplots(figsize=(5, 2.6))
                    ax.plot(sweep["k"], sweep[metric], marker="o", color="#1f77b4")
                    ax.set_title(title, fontsize=10)
                    ax.set_xlabel("Number of groups")
                    show(fig)
            best = sweep.loc[sweep["silhouette"].idxmax()]
            st.caption(f"The best silhouette is at {int(best['k'])} groups ({best['silhouette']:.3f}).")

        st.subheader("Speed: all genes vs PCA first")
        if st.button("Run speed comparison"):
            with st.spinner("Timing each step both ways..."):
                s = result["settings"]
                st.session_state.speed = pl.speed_comparison(result["Z"], s["n_pca"], s["k"],
                                                             s["method"], s["perplexity"])
        if "speed" in st.session_state:
            table, ari = st.session_state.speed
            left, right = st.columns([3, 2])
            with left:
                st.dataframe(table.round(3).rename(columns={
                    "step": "Step", "all_genes_seconds": "All genes (s)", "pca_first_seconds": "PCA first (s)",
                    "all_genes_silhouette": "Silhouette, all genes", "pca_first_silhouette": "Silhouette, PCA first",
                    "speedup": "Speedup"}), hide_index=True)
                st.caption(f"PCA times include fitting PCA itself. The two groupings agree with an "
                           f"adjusted Rand index of {ari:.2f} (1 means identical).")
            with right:
                fig, ax = plt.subplots(figsize=(4.5, 2.8))
                x = np.arange(len(table))
                ax.bar(x - 0.2, table["all_genes_seconds"], 0.4, label="All genes", color="#7f7f7f")
                ax.bar(x + 0.2, table["pca_first_seconds"], 0.4, label="PCA first", color="#1f77b4")
                ax.set_xticks(x, table["step"], fontsize=8)
                ax.set_ylabel("Seconds")
                ax.legend(fontsize=8, frameon=False)
                show(fig)

        st.subheader("Genes that define each group")
        markers = st.session_state.classifier.markers
        g = st.selectbox("Group", range(len(markers)), format_func=pl.group_name)
        st.dataframe(markers[g].head(10).rename(columns={"difference": "difference (std devs above the rest)"}),
                     hide_index=True)

        assignments = pd.DataFrame({
            "patient_id": result["patient_ids"],
            "group": [pl.group_name(g) for g in labels],
            "pc1": result["pca_coords"][:, 0], "pc2": result["pca_coords"][:, 1],
            "tsne1": result["tsne_coords"][:, 0], "tsne2": result["tsne_coords"][:, 1],
        })
        summary = pd.Series({**result["settings"], **m, **{f"{key}_seconds": v for key, v in t.items()}})
        left, right = st.columns(2)
        left.download_button("Download group assignments (CSV)", assignments.to_csv(index=False),
                             "group_assignments.csv", "text/csv")
        right.download_button("Download metrics (CSV)", summary.rename_axis("name").to_frame("value").to_csv(),
                              "metrics.csv", "text/csv")
    st.caption(FOOTER)

# Tab 2: unseen patients

with tab_classify:
    if result is None:
        st.info("Run the pipeline first. The groups it finds are what new patients are matched against.")
    else:
        clf, X_all = st.session_state.classifier, st.session_state.X
        source = st.radio("Patient", ["Unseen patient from the dataset", "Import a patient file"], horizontal=True)
        pid, profile = None, None
        if source == "Import a patient file":
            upload = st.file_uploader("CSV with gene names in the first column and one column per patient",
                                      type=["csv"], key="patient_file")
            if upload is not None:
                try:
                    imported = load_upload(upload.getvalue())
                except ValueError as e:
                    st.error(str(e))
                else:
                    pid = st.selectbox(f"Imported patient ({len(imported)} in the file)", list(imported.index))
                    profile = imported.loc[pid]
        elif result["holdout_ids"]:
            st.write(f"{len(result['holdout_ids'])} patients from the dataset were kept out of the grouping. "
                     "Pick one to see which group the model puts them in.")
            pid = st.selectbox("Unseen patient", result["holdout_ids"])
            profile = X_all.loc[pid]
        else:
            st.info("No patients were held back. Raise Patients held back for testing in the sidebar and run again.")
        if st.button("Classify", type="primary", disabled=profile is None):
            try:
                st.session_state.classified = (pid, clf.classify(profile))
            except ValueError as e:
                st.error(str(e))

        if "classified" in st.session_state:
            shown_id, res = st.session_state.classified
            g = res["group"]
            left, right = st.columns([2, 3])
            with left:
                st.caption(f"Result for {shown_id}")
                st.markdown(f"## {pl.group_name(g)}")
                st.write(f"**{res['confidence_label']}**, {res['confidence']:.0%} of the nearest-patient vote.")
                if res["n_missing_genes"]:
                    st.warning(f"{res['n_missing_genes']} genes were missing and filled with the dataset median.")
                if res["n_extra_genes"]:
                    st.warning(f"{res['n_extra_genes']} genes are not in the training data and were ignored.")
                proba = res["probabilities"]
                fig, ax = plt.subplots(figsize=(4.5, 0.45 * len(proba) + 0.6))
                names = [pl.group_name(i) for i in range(len(proba))]
                ax.barh(names, proba, color=pl.GROUP_COLORS[:len(proba)])
                ax.invert_yaxis()
                ax.set_xlim(0, 1)
                ax.set_xlabel("Share of the vote")
                for i, p in enumerate(proba):
                    ax.text(min(p + 0.02, 0.9), i, f"{p:.0%}", va="center", fontsize=8)
                show(fig)
            with right:
                fig, ax = plt.subplots(figsize=(5.5, 4.2))
                scatter(ax, result["tsne_coords"], result["labels"], alpha=0.3)
                x, y = res["tsne_xy"]
                ax.scatter([x], [y], marker="*", s=380, color="black", edgecolors="white", linewidths=1.2,
                           zorder=5)
                ax.annotate(shown_id, (x, y), xytext=(8, 8), textcoords="offset points", fontsize=8)
                ax.set_title("Where this patient lands on the t-SNE map", fontsize=11)
                ax.set_xticks([])
                ax.set_yticks([])
                show(fig)

            left, right = st.columns(2)
            with left:
                st.subheader("Most similar patients")
                st.dataframe(res["neighbors"], hide_index=True)
            with right:
                st.subheader(f"Genes behind {pl.group_name(g)}")
                st.dataframe(res["driving_genes"], hide_index=True)
                st.caption(f"Genes that mark {pl.group_name(g)}, with this patient's level "
                           "in standard deviations from the average patient.")

        st.subheader("All unseen patients")
        if st.button("Classify all of them", disabled=not result["holdout_ids"]):
            rows = []
            for p in result["holdout_ids"]:
                r = clf.classify(X_all.loc[p])
                rows.append({"patient_id": p, "group": pl.group_name(r["group"]),
                             "confidence": round(r["confidence"], 2), "match": r["confidence_label"]})
            st.session_state.holdout_table = pd.DataFrame(rows)
        if "holdout_table" in st.session_state:
            table = st.session_state.holdout_table
            counts = table["match"].value_counts()
            st.caption(", ".join(f"{counts.get(lbl, 0)} {lbl.lower()}"
                                 for lbl in ("Clear match", "Moderate match", "Borderline")) + ".")
            st.dataframe(table, hide_index=True)
            st.download_button("Download predictions (CSV)", table.to_csv(index=False),
                               "unseen_patient_groups.csv", "text/csv")
    st.caption(FOOTER)
