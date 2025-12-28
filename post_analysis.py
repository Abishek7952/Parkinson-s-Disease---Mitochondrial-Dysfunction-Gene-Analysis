import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import GEOparse
import os

# ==============================
# SETUP
# ==============================
os.makedirs("outputs", exist_ok=True)
sns.set(style="whitegrid")

print("Loading datasets...")

# ------------------------------
# Load expression data
# ------------------------------
expression_df = pd.read_csv(
    "GSE6613_series_matrix .txt",
    sep="\t",
    comment="!",
    header=0
)

# ------------------------------
# Load MitoCarta
# ------------------------------
mitocarta_df = pd.read_excel(
    "mitocartagenes.xls",
    sheet_name=1
)
mito_genes = mitocarta_df["Symbol"].dropna().unique()

# ------------------------------
# Load GPL96 annotation
# ------------------------------
gpl = GEOparse.get_GEO("GPL96", destdir=".")
annotation_df = gpl.table[["ID", "Gene Symbol"]]
annotation_df.columns = ["ID_REF", "GeneSymbol"]
annotation_df = annotation_df.dropna()
annotation_df = annotation_df[annotation_df["GeneSymbol"] != ""]

# ------------------------------
# Merge expression + annotation
# ------------------------------
expr_annotated = expression_df.merge(
    annotation_df,
    on="ID_REF",
    how="left"
)

mito_expr = expr_annotated[
    expr_annotated["GeneSymbol"].isin(mito_genes)
].copy()

# ------------------------------
# Define sample groups
# ------------------------------
sample_cols = expression_df.columns[1:]
control_samples = sample_cols[:50]
pd_samples = sample_cols[50:]

# ------------------------------
# Calculate means & fold change
# ------------------------------
mito_expr["Control_Mean"] = mito_expr[control_samples].mean(axis=1)
mito_expr["PD_Mean"] = mito_expr[pd_samples].mean(axis=1)
mito_expr["FoldChange_PD_vs_Control"] = mito_expr["PD_Mean"] / mito_expr["Control_Mean"]

mito_expr.replace([np.inf, -np.inf], np.nan, inplace=True)
mito_expr.dropna(subset=["FoldChange_PD_vs_Control"], inplace=True)

# ==============================
# EXPORT CSV FILES
# ==============================
print("Exporting CSV files...")

ranked_genes = mito_expr.sort_values(
    by="FoldChange_PD_vs_Control",
    ascending=False
)

ranked_genes.to_csv("outputs/ranked_mito_genes.csv", index=False)

# Gene-level averaged matrix
gene_level_expr = mito_expr.groupby("GeneSymbol")[control_samples.tolist() + pd_samples.tolist()].mean()
gene_level_expr.to_csv("outputs/mito_expression_matrix.csv")

# ==============================
# HEATMAP (TOP 20 GENES)
# ==============================
print("Generating heatmap...")

top20_genes = ranked_genes["GeneSymbol"].unique()[:20]
heatmap_data = gene_level_expr.loc[top20_genes]

plt.figure(figsize=(14, 8))
sns.heatmap(
    heatmap_data,
    cmap="coolwarm",
    yticklabels=True
)
plt.title("Top 20 Differentially Expressed Mitochondrial Genes")
plt.xlabel("Samples")
plt.ylabel("Genes")
plt.tight_layout()
plt.savefig("outputs/heatmap_top20.png", dpi=300)
plt.close()

# ==============================
# VOLCANO PLOT
# ==============================
print("Generating volcano plot...")

gene_fc = ranked_genes.groupby("GeneSymbol")["FoldChange_PD_vs_Control"].mean()
log_fc = np.log2(gene_fc)

mean_expr = np.log10(
    gene_level_expr[pd_samples].mean(axis=1) + 1
)

plt.figure(figsize=(8, 6))
plt.scatter(mean_expr, log_fc, alpha=0.6)
plt.axhline(1, linestyle="--", color="grey")
plt.axhline(-1, linestyle="--", color="grey")
plt.xlabel("log10 Mean Expression (PD)")
plt.ylabel("log2 Fold Change (PD / Control)")
plt.title("Volcano Plot of Mitochondrial Genes")
plt.tight_layout()
plt.savefig("outputs/volcano_plot.png", dpi=300)
plt.close()

# ==============================
# BOXPLOT (TOP GENE – GENE LEVEL)
# ==============================
print("Generating boxplot...")

top_gene = top20_genes[0]
top_gene_expr = gene_level_expr.loc[top_gene]

box_data = pd.DataFrame({
    "Expression": list(top_gene_expr[control_samples]) +
                  list(top_gene_expr[pd_samples]),
    "Group": ["Control"] * len(control_samples) +
             ["Parkinson's"] * len(pd_samples)
})

plt.figure(figsize=(6, 5))
sns.boxplot(x="Group", y="Expression", data=box_data)
plt.title(f"Expression of {top_gene}")
plt.tight_layout()
plt.savefig("outputs/boxplot_top_gene.png", dpi=300)
plt.close()

print("\nPOST-ANALYSIS COMPLETE")
print("All files saved in /outputs/")
