import pandas as pd
import numpy as np
import GEOparse

print("\n==============================")
print("STEP 1: LOAD GSE6613 DATASET")
print("==============================")

# -----------------------------
# Load GSE6613 expression data
# -----------------------------
expression_file = "GSE6613_series_matrix .txt"

expression_df = pd.read_csv(
    expression_file,
    sep="\t",
    comment="!",
    header=0
)

print("Expression dataset shape:")
print(expression_df.shape)


print("\n==============================")
print("STEP 2: LOAD MITOCARTA GENE LIST")
print("==============================")

# -----------------------------
# Load MitoCarta mitochondrial gene list
# -----------------------------
mitocarta_df = pd.read_excel(
    "mitocartagenes.xls",
    sheet_name=1   # A Human MitoCarta3.0
)

mito_genes = mitocarta_df["Symbol"].dropna().unique()
print(f"Total mitochondrial genes: {len(mito_genes)}")


print("\n==============================")
print("STEP 3: LOAD GPL96 PROBE ANNOTATION")
print("==============================")

# -----------------------------
# Load platform annotation using GEOparse
# -----------------------------
gpl = GEOparse.get_GEO("GPL96", destdir=".")

annotation_df = gpl.table[["ID", "Gene Symbol"]]
annotation_df.columns = ["ID_REF", "GeneSymbol"]

annotation_df = annotation_df.dropna()
annotation_df = annotation_df[annotation_df["GeneSymbol"] != ""]

print("Annotation dataset shape:")
print(annotation_df.shape)


print("\n==============================")
print("STEP 4: MAP PROBES TO GENE SYMBOLS")
print("==============================")

# Merge expression data with annotation
expr_annotated = expression_df.merge(
    annotation_df,
    on="ID_REF",
    how="left"
)

print("Annotated expression shape:")
print(expr_annotated.shape)


print("\n==============================")
print("STEP 5: FILTER MITOCHONDRIAL GENES")
print("==============================")

mito_expr = expr_annotated[
    expr_annotated["GeneSymbol"].isin(mito_genes)
]

print("Mitochondrial expression shape:")
print(mito_expr.shape)


print("\n==============================")
print("STEP 6: DEFINE PD VS CONTROL GROUPS")
print("==============================")

# According to GSE6613 metadata
sample_cols = mito_expr.columns[1:-1]  # exclude ID_REF & GeneSymbol

control_samples = sample_cols[:50]
pd_samples = sample_cols[50:]

print(f"Control samples: {len(control_samples)}")
print(f"Parkinson's samples: {len(pd_samples)}")


print("\n==============================")
print("STEP 7: CALCULATE MEAN EXPRESSION & FOLD CHANGE")
print("==============================")

mito_expr["Control_Mean"] = mito_expr[control_samples].mean(axis=1)
mito_expr["PD_Mean"] = mito_expr[pd_samples].mean(axis=1)

mito_expr["FoldChange_PD_vs_Control"] = (
    mito_expr["PD_Mean"] / mito_expr["Control_Mean"]
)

mito_expr = mito_expr.replace([np.inf, -np.inf], np.nan)
mito_expr = mito_expr.dropna(subset=["FoldChange_PD_vs_Control"])


print("\n==============================")
print("STEP 8: RANK DIFFERENTIALLY EXPRESSED MITO GENES")
print("==============================")

ranked = mito_expr.sort_values(
    by="FoldChange_PD_vs_Control",
    ascending=False
)

print("\nTop 10 UPREGULATED mitochondrial genes (PD):")
print(ranked[["GeneSymbol", "FoldChange_PD_vs_Control"]].head(10))

print("\nTop 10 DOWNREGULATED mitochondrial genes (PD):")
print(ranked[["GeneSymbol", "FoldChange_PD_vs_Control"]].tail(10))


print("\n==============================")
print("ANALYSIS COMPLETE")
print("==============================")
