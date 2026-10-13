from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

PROJECT_DIR = Path(__file__).resolve().parents[1]

targetscan_path = (
    PROJECT_DIR
    / "data"
    / "raw"
    / "TargetScan8.0__miR-186-5p.Human.predicted_targets.txt"
)

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

print("1. Загружаем предсказанные мишени TargetScan...")

targets = pd.read_csv(
    targetscan_path,
    sep="\t",
    low_memory=False
)

print("Размер таблицы TargetScan:", targets.shape)
print("Названия первых столбцов:")
print(targets.columns[:5].tolist())

target_column = targets.columns[0]

target_genes = (
    targets[target_column]
    .astype(str)
    .str.strip()
    .replace("nan", np.nan)
    .dropna()
    .unique()
    .tolist()
)

print("Число уникальных генов-мишеней:", len(target_genes))
print("Первые 10 мишеней:", target_genes[:10])

print("\n2. Загружаем матрицы экспрессии...")

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

expression = pd.read_csv(
    expression_path,
    sep="\t",
    compression="gzip"
)

expression = expression.rename(
    columns={"sample": "gene_name"}
)

expression["gene_name"] = (
    expression["gene_name"]
    .astype(str)
    .str.strip()
)

expression = expression.drop_duplicates(
    subset="gene_name"
).set_index("gene_name")

expression = expression[common_samples]

available_targets = [
    gene
    for gene in target_genes
    if gene in expression.index
]

print("Мишеней, найденных в матрице TCGA-SKCM:",
      len(available_targets))

anxa2p1 = pseudogene_expression.loc[
    "ANXA2P1",
    common_samples
]

mir186 = mirna_expression.loc[
    "MIMAT0000456",
    common_samples
]

results = []

print("\n3. Считаем корреляции...")

for gene in available_targets:
    target_expression = expression.loc[
        gene,
        common_samples
    ]

    r_mirna_target, p_mirna_target = spearmanr(
        mir186,
        target_expression
    )

    r_pseudogene_target, p_pseudogene_target = spearmanr(
        anxa2p1,
        target_expression
    )

    results.append({
        "pseudogene": "ANXA2P1",
        "mirna": "hsa-miR-186-5p",
        "target_gene": gene,
        "r_mirna_target": r_mirna_target,
        "p_mirna_target": p_mirna_target,
        "r_pseudogene_target": r_pseudogene_target,
        "p_pseudogene_target": p_pseudogene_target
    })

results = pd.DataFrame(results)

results["FDR_mirna_target"] = np.nan
results["FDR_pseudogene_target"] = np.nan

valid_mirna = results["p_mirna_target"].notna()
valid_pseudogene = results["p_pseudogene_target"].notna()

results.loc[
    valid_mirna,
    "FDR_mirna_target"
] = multipletests(
    results.loc[valid_mirna, "p_mirna_target"],
    method="fdr_bh"
)[1]

results.loc[
    valid_pseudogene,
    "FDR_pseudogene_target"
] = multipletests(
    results.loc[valid_pseudogene, "p_pseudogene_target"],
    method="fdr_bh"
)[1]

results = results.sort_values(
    ["r_mirna_target", "r_pseudogene_target"],
    ascending=[True, False]
)

all_results_path = (
    tables_dir
    / "anxa2p1_mir186_all_targetscan_targets_correlations.tsv"
)

results.to_csv(
    all_results_path,
    sep="\t",
    index=False
)

cerna_candidates = results[
    (results["r_mirna_target"] <= -0.15)
    & (results["FDR_mirna_target"] < 0.05)
    & (results["r_pseudogene_target"] >= 0.15)
    & (results["FDR_pseudogene_target"] < 0.05)
].copy()

cerna_candidates = cerna_candidates.sort_values(
    ["r_mirna_target", "r_pseudogene_target"],
    ascending=[True, False]
)

candidate_path = (
    tables_dir
    / "anxa2p1_mir186_potential_cerna_targets.tsv"
)

cerna_candidates.to_csv(
    candidate_path,
    sep="\t",
    index=False
)

print("\nВсего проверено мишеней:", len(results))
print("Потенциальных ceRNA-троек:", len(cerna_candidates))

print("\nКандидаты:")
if len(cerna_candidates) == 0:
    print("Строгим критериям ни одна мишень не соответствует.")
else:
    print(
        cerna_candidates[
            [
                "target_gene",
                "r_mirna_target",
                "FDR_mirna_target",
                "r_pseudogene_target",
                "FDR_pseudogene_target"
            ]
        ].head(20).to_string(index=False)
    )

print("\nСохранены файлы:")
print(all_results_path)
print(candidate_path)