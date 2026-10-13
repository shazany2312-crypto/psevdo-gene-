from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

PROJECT_DIR = Path(__file__).resolve().parents[1]

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

print("Псевдогенов:", pseudogene_expression.shape[0])
print("МикроРНК:", mirna_expression.shape[0])
print("Общих образцов:", len(common_samples))

candidates = ["ANXA2P1", "ANXA2P2", "ANXA2P3"]

results = []

for pseudogene in candidates:
    pseudo_values = pseudogene_expression.loc[
        pseudogene,
        common_samples
    ]

    for mirna_id in mirna_expression.index:
        mirna_values = mirna_expression.loc[
            mirna_id,
            common_samples
        ]

        correlation, pvalue = spearmanr(
            pseudo_values,
            mirna_values
        )

        results.append({
            "pseudogene": pseudogene,
            "mirna_id": mirna_id,
            "spearman_r": correlation,
            "pvalue": pvalue
        })

results = pd.DataFrame(results)

print("\nПар с неопределённой корреляцией или p-value:",
      results["pvalue"].isna().sum())

results["FDR"] = float("nan")

valid_pvalues = results["pvalue"].notna()

results.loc[valid_pvalues, "FDR"] = multipletests(
    results.loc[valid_pvalues, "pvalue"],
    method="fdr_bh"
)[1]

results = results.sort_values(
    ["pseudogene", "spearman_r"],
    ascending=[True, True]
)

results = results.sort_values(
    ["pseudogene", "spearman_r"],
    ascending=[True, True]
)

all_results_path = (
    tables_dir
    / "anxa2_pseudogene_mirna_correlations.tsv"
)

results.to_csv(
    all_results_path,
    sep="\t",
    index=False
)

negative_candidates = results[
    (results["spearman_r"] <= -0.20)
    & (results["FDR"] < 0.05)
].copy()

negative_candidates_path = (
    tables_dir
    / "anxa2_pseudogene_negative_mirna_candidates.tsv"
)

negative_candidates.to_csv(
    negative_candidates_path,
    sep="\t",
    index=False
)

print("\nВсего пар псевдоген–микроРНК:", len(results))
print("Кандидатов с r ≤ -0.20 и FDR < 0.05:",
      len(negative_candidates))

for pseudogene in candidates:
    print(f"\nТоп-10 наиболее отрицательных корреляций: {pseudogene}")

    top = results[
        results["pseudogene"] == pseudogene
    ].head(10)

    print(
        top[
            [
                "mirna_id",
                "spearman_r",
                "pvalue",
                "FDR"
            ]
        ].to_string(index=False)
    )

print("\nСохранены таблицы:")
print(all_results_path)
print(negative_candidates_path)