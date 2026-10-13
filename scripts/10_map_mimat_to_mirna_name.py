from pathlib import Path
import re
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]

fasta_path = (
    PROJECT_DIR
    / "data"
    / "raw"
    / "mature.fa"
)

correlations_path = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "anxa2_pseudogene_mirna_correlations.tsv"
)

tables_dir = PROJECT_DIR / "results" / "tables"

mirna_map = []

with open(fasta_path, "r", encoding="utf-8") as fasta_file:
    for line in fasta_file:
        if not line.startswith(">"):
            continue

        parts = line.strip()[1:].split()

        mirna_name = parts[0]

        mirna_id = next(
            (
                part
                for part in parts
                if part.startswith("MIMAT")
            ),
            None
        )

        if mirna_id is not None:
            mirna_map.append({
                "mirna_name": mirna_name,
                "mirna_id": mirna_id
            })

mirna_map = pd.DataFrame(
    mirna_map,
    columns=["mirna_name", "mirna_id"]
)

human_mirna_map = mirna_map[
    mirna_map["mirna_name"].str.startswith("hsa-")
].copy()

human_mirna_map = human_mirna_map.drop_duplicates(
    subset="mirna_id"
)

print("Всего записей зрелых микроРНК:", len(mirna_map))
print("Записей микроРНК человека:", len(human_mirna_map))

print("\nПример соответствий:")
print(human_mirna_map.head(10).to_string(index=False))

correlations = pd.read_csv(
    correlations_path,
    sep="\t"
)

correlations_named = correlations.merge(
    human_mirna_map,
    on="mirna_id",
    how="left"
)

columns = [
    "pseudogene",
    "mirna_id",
    "mirna_name",
    "spearman_r",
    "pvalue",
    "FDR"
]

correlations_named = correlations_named[columns]

output_path = (
    tables_dir
    / "anxa2_pseudogene_mirna_correlations_named.tsv"
)

correlations_named.to_csv(
    output_path,
    sep="\t",
    index=False
)

negative_candidates = correlations_named[
    (correlations_named["spearman_r"] <= -0.20)
    & (correlations_named["FDR"] < 0.05)
].copy()

candidate_path = (
    tables_dir
    / "anxa2_pseudogene_negative_mirna_candidates_named.tsv"
)

negative_candidates.to_csv(
    candidate_path,
    sep="\t",
    index=False
)

print("\nКандидаты с понятными названиями:")
print(
    negative_candidates[
        [
            "pseudogene",
            "mirna_name",
            "mirna_id",
            "spearman_r",
            "FDR"
        ]
    ].to_string(index=False)
)

print("\nСохранены файлы:")
print(output_path)
print(candidate_path)