from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr

PROJECT_DIR = Path(__file__).resolve().parents[1]

expression_path = (
    PROJECT_DIR
    / "data"
    / "raw"
    / "TCGA_SKCM_sampleMap_HiSeqV2.gz"
)

pseudogene_path = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "skcm_pseudogene_common_samples.tsv.gz"
)

mirna_path = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "skcm_mirna_common_samples.tsv.gz"
)

tables_dir = PROJECT_DIR / "results" / "tables"

pseudogene_expression = pd.read_csv(
    pseudogene_path,
    sep="\t",
    compression="gzip",
    index_col=0
)

mirna_expression = pd.read_csv(
    mirna_path,
    sep="\t",
    compression="gzip",
    index_col=0
)

common_samples = pseudogene_expression.columns.tolist()

mirna_expression = mirna_expression[common_samples]

print("Общих образцов:", len(common_samples))

target_genes = ["CYLD", "PTPRC", "IL7R"]

expression = pd.read_csv(
    expression_path,
    sep="\t",
    compression="gzip"
)

expression = expression.rename(
    columns={"sample": "gene_name"}
)

expression = expression.set_index("gene_name")

available_targets = [
    gene
    for gene in target_genes
    if gene in expression.index
]

print("Кандидатных мРНК в матрице:", available_targets)

expression = expression[common_samples]

anxa2p1 = pseudogene_expression.loc[
    "ANXA2P1",
    common_samples
]

mir186 = mirna_expression.loc[
    "MIMAT0000456",
    common_samples
]

results = []

for gene in available_targets:
    target_expression = expression.loc[
        gene,
        common_samples
    ]

    r_mirna_target, p_mirna_target = spearmanr(
        mir186,
        target_expression
    )

    r_pseudo_target, p_pseudo_target = spearmanr(
        anxa2p1,
        target_expression
    )

    results.append({
        "pseudogene": "ANXA2P1",
        "mirna": "hsa-miR-186-5p",
        "target_gene": gene,
        "r_mirna_target": r_mirna_target,
        "p_mirna_target": p_mirna_target,
        "r_pseudogene_target": r_pseudo_target,
        "p_pseudogene_target": p_pseudo_target
    })

results = pd.DataFrame(results)

output_path = (
    tables_dir
    / "anxa2p1_mir186_candidate_target_correlations.tsv"
)

results.to_csv(
    output_path,
    sep="\t",
    index=False
)

print("\nКорреляции потенциальных ceRNA-троек:")
print(results.to_string(index=False))

print("\nСохранено:")
print(output_path)