from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

PROJECT_DIR = Path(__file__).resolve().parents[1]

expression_path = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "skcm_tumor_expressed_processed_pseudogenes.tsv.gz"
)

annotation_path = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "skcm_sample_annotation.tsv"
)

tables_dir = PROJECT_DIR / "results" / "tables"

print("1. Загружаем экспрессию псевдогенов...")
expression = pd.read_csv(
    expression_path,
    sep="\t",
    compression="gzip",
    index_col=0
)

print("2. Загружаем аннотацию образцов...")
annotation = pd.read_csv(
    annotation_path,
    sep="\t"
)

primary_samples = annotation.loc[
    annotation["sample_type"] == "Primary tumor",
    "sample_id"
].tolist()

metastatic_samples = annotation.loc[
    annotation["sample_type"].isin([
        "Metastatic tumor",
        "Additional metastatic tumor"
    ]),
    "sample_id"
].tolist()

print("Первичных опухолей:", len(primary_samples))
print("Метастатических опухолей:", len(metastatic_samples))

results = []

for pseudogene in expression.index:
    primary_values = expression.loc[
        pseudogene,
        primary_samples
    ].dropna()

    metastatic_values = expression.loc[
        pseudogene,
        metastatic_samples
    ].dropna()

    statistic, pvalue = mannwhitneyu(
        primary_values,
        metastatic_values,
        alternative="two-sided"
    )

    log2fc = (
        metastatic_values.mean()
        - primary_values.mean()
    )

    results.append({
        "pseudogene": pseudogene,
        "mean_primary": primary_values.mean(),
        "mean_metastatic": metastatic_values.mean(),
        "median_primary": primary_values.median(),
        "median_metastatic": metastatic_values.median(),
        "log2FC_metastatic_vs_primary": log2fc,
        "pvalue": pvalue
    })

results = pd.DataFrame(results)

results["FDR"] = multipletests(
    results["pvalue"],
    method="fdr_bh"
)[1]

results = results.sort_values(
    "FDR",
    ascending=True
)

significant = results[
    (results["FDR"] < 0.05)
    & (results["log2FC_metastatic_vs_primary"].abs() > 0.5)
].copy()

all_results_path = (
    tables_dir
    / "primary_vs_metastatic_pseudogene_expression.tsv"
)

significant_path = (
    tables_dir
    / "significant_pseudogenes_primary_vs_metastatic.tsv"
)

results.to_csv(
    all_results_path,
    sep="\t",
    index=False
)

significant.to_csv(
    significant_path,
    sep="\t",
    index=False
)

print("\nВсего проверено псевдогенов:", len(results))
print("Статистически значимых кандидатов:", len(significant))

print("\nНаиболее значимые различия:")
print(
    results[
        [
            "pseudogene",
            "mean_primary",
            "mean_metastatic",
            "log2FC_metastatic_vs_primary",
            "pvalue",
            "FDR"
        ]
    ].head(10).to_string(index=False)
)

print("\nСохранены файлы:")
print(all_results_path)
print(significant_path)