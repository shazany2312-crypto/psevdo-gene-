from pathlib import Path
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]

expression_path = (
    PROJECT_DIR
    / "data"
    / "raw"
    / "TCGA_SKCM_sampleMap_HiSeqV2.gz"
)

pseudogene_path = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "processed_pseudogenes.tsv"
)

processed_dir = PROJECT_DIR / "data" / "processed"
processed_dir.mkdir(parents=True, exist_ok=True)

print("1. Загружаем матрицу экспрессии меланомы...")
expression = pd.read_csv(
    expression_path,
    sep="\t",
    compression="gzip"
)

print("Размер исходной матрицы:", expression.shape)

expression = expression.rename(
    columns={"sample": "gene_name"}
)

expression["gene_name"] = expression["gene_name"].astype(str).str.strip()

print("2. Загружаем список процессированных псевдогенов...")
pseudogenes = pd.read_csv(
    pseudogene_path,
    sep="\t"
)

pseudogenes["gene_name"] = (
    pseudogenes["gene_name"]
    .astype(str)
    .str.strip()
)

pseudogene_names = set(
    pseudogenes["gene_name"]
    .dropna()
)

print("Всего псевдогенов в аннотации:", len(pseudogene_names))

pseudogene_expression = expression[
    expression["gene_name"].isin(pseudogene_names)
].copy()

pseudogene_expression = pseudogene_expression.drop_duplicates(
    subset="gene_name"
)

pseudogene_expression = pseudogene_expression.set_index(
    "gene_name"
)

print("Псевдогенов найдено в RNA-seq-матрице:",
      pseudogene_expression.shape[0])

print("Число образцов:", pseudogene_expression.shape[1])

output_path = (
    processed_dir
    / "skcm_processed_pseudogene_expression.tsv.gz"
)

pseudogene_expression.to_csv(
    output_path,
    sep="\t",
    compression="gzip"
)

overlap_names = pseudogene_expression.index.tolist()

overlap_table = pseudogenes[
    pseudogenes["gene_name"].isin(overlap_names)
].copy()

overlap_table_path = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "pseudogenes_found_in_skcm.tsv"
)

overlap_table.to_csv(
    overlap_table_path,
    sep="\t",
    index=False
)

print("\nРезультаты сохранены:")
print(output_path)
print(overlap_table_path)

print("\nПервые 10 найденных псевдогенов:")
print(pseudogene_expression.index[:10].tolist())