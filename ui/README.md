# Genomic Subtype Clustering app

A Streamlit app for the project pipeline: standardize, PCA, t-SNE, then K-Means or hierarchical clustering.

It has two tabs.

**Run pipeline.** Pick the PCA components, t-SNE perplexity, grouping method, number of groups and the space to cluster on, then run. You get silhouette, Davies-Bouldin and Calinski-Harabasz scores, the PCA and t-SNE maps, explained variance, a comparison of 2 to 8 groups, a timing comparison of all genes against PCA first, the genes that define each group, and CSV downloads.

**Classify a patient.** A share of the TCGA patients (20% by default) is held back before grouping. Pick one of these unseen patients and the app places them in a group using the 15 most similar patients in PCA space. It shows the confidence, where they land on the t-SNE map, their closest patients and the genes behind the match. You can also classify every unseen patient at once.

## Run

From the `ML_PROJECT` folder:

```
pip install -r ui/requirements.txt
streamlit run ui/app.py
```

The built-in data is `RNA_seq_1.csv` (TCGA, 528 patients, 4,571 genes). To use other data, upload a CSV in the same layout: the first column holds gene names, and every other column is one patient.
