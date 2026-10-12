# Test data

All files are cut from the real TCGA data in `RNA_seq_1.csv`, and nothing in them is made up. They use the same layout as the original: genes as rows, one column per patient.

- `test_set.csv`: the 119 patients the apps hold back (20% split, seed 42), with all 4,571 genes. The model never sees these patients when it builds the groups.
- `small_sample.csv`: 60 patients and the 300 most variable genes (0.13 MB). In the Streamlit app, upload it under "Upload my own CSV" for a pipeline run that takes under a second.
- `extra_patient.csv`: one held-back patient on its own (gene, value columns). Import it with "Import patient file" in the galaxy or "Import a patient file" in the Streamlit Classify tab.

To rebuild the files: `python test_data/make_test_data.py`
