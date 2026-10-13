from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

PROJECT_DIR = Path(__file__).resolve().parents[1]

results_path = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "primary_vs_metastatic_pseudogene_expression.tsv"
)

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

figures_dir = PROJECT_DIR / "figures"
figures_dir.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="notebook")

results = pd.read_csv(
    results_path,
    sep="\t"
)

results["minus_log10_FDR"] = -np.log10(
    results["FDR"].clip(lower=1e-300)
)

results["significant"] = (
    (results["FDR"] < 0.05)
    & (results["log2FC_metastatic_vs_primary"].abs() > 0.5)
)

plt.figure(figsize=(9, 6))

sns.scatterplot(
    data=results,
    x="log2FC_metastatic_vs_primary",
    y="minus_log10_FDR",
    hue="significant",
    palette={False: "gray", True: "crimson"},
    s=70,
    legend=False
)

plt.axvline(-0.5, color="black", linestyle="--", linewidth=1)
plt.axvline(0.5, color="black", linestyle="--", linewidth=1)
plt.axhline(
    -np.log10(0.05),
    color="black",
    linestyle="--",
    linewidth=1
)

for _, row in results[results["significant"]].iterrows():
    plt.text(
        row["log2FC_metastatic_vs_primary"] + 0.02,
        row["minus_log10_FDR"] + 0.03,
        row["pseudogene"],
        fontsize=9
    )

plt.xlabel(
    "Разность средней экспрессии: метастазы − первичные опухоли"
)
plt.ylabel("−log10(FDR)")
plt.title(
    "Различия экспрессии процессированных псевдогенов\n"
    "между первичной и метастатической меланомой"
)

plt.tight_layout()

volcano_path = figures_dir / "volcano_primary_vs_metastatic.png"

plt.savefig(
    volcano_path,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

expression = pd.read_csv(
    expression_path,
    sep="\t",
    compression="gzip",
    index_col=0
)

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

candidates = ["ANXA2P1", "ANXA2P2", "ANXA2P3"]

plot_data = []

for pseudogene in candidates:
    for sample in primary_samples:
        plot_data.append({
            "pseudogene": pseudogene,
            "sample_type": "Первичная опухоль",
            "expression": expression.loc[pseudogene, sample]
        })

    for sample in metastatic_samples:
        plot_data.append({
            "pseudogene": pseudogene,
            "sample_type": "Метастаз",
            "expression": expression.loc[pseudogene, sample]
        })

plot_data = pd.DataFrame(plot_data)

plt.figure(figsize=(11, 6))

sns.boxplot(
    data=plot_data,
    x="pseudogene",
    y="expression",
    hue="sample_type",
    palette={
        "Первичная опухоль": "#4C72B0",
        "Метастаз": "#DD8452"
    },
    showfliers=False
)

sns.stripplot(
    data=plot_data,
    x="pseudogene",
    y="expression",
    hue="sample_type",
    dodge=True,
    color="black",
    alpha=0.12,
    size=2,
    legend=False
)

plt.xlabel("Процессированный псевдоген")
plt.ylabel("log2(norm_count + 1)")
plt.title(
    "Экспрессия значимых процессированных псевдогенов\n"
    "в первичной и метастатической меланоме"
)

plt.legend(
    title="Тип образца",
    loc="best"
)

plt.tight_layout()

boxplot_path = figures_dir / "boxplots_anxa2_pseudogenes.png"

plt.savefig(
    boxplot_path,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Графики сохранены:")
print(volcano_path)
print(boxplot_path)