from pathlib import Path
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]

input_path = (
    PROJECT_DIR
    / "data"
    / "processed"
    / "skcm_processed_pseudogene_expression.tsv.gz"
)

processed_dir = PROJECT_DIR / "data" / "processed"
tables_dir = PROJECT_DIR / "results" / "tables"

print("1. Загружаем матрицу псевдогенов...")
expression = pd.read_csv(
    input_path,
    sep="\t",
    compression="gzip",
    index_col=0
)

print("Размер матрицы:", expression.shape)

sample_ids = expression.columns.tolist()

sample_annotation = pd.DataFrame({
    "sample_id": sample_ids
})

sample_annotation["sample_type_code"] = (
    sample_annotation["sample_id"]
    .str.split("-")
    .str[-1]
    .str[:2]
)

sample_type_names = {
    "01": "Primary tumor",
    "06": "Metastatic tumor",
    "07": "Additional metastatic tumor",
    "11": "Solid tissue normal"
}

sample_annotation["sample_type"] = (
    sample_annotation["sample_type_code"]
    .map(sample_type_names)
    .fillna("Other")
)

print("\nКоличество образцов по типам:")
print(sample_annotation["sample_type"].value_counts())

sample_annotation_path = (
    tables_dir
    / "skcm_sample_annotation.tsv"
)

sample_annotation.to_csv(
    sample_annotation_path,
    sep="\t",
    index=False
)

tumor_codes = ["01", "06", "07"]

tumor_samples = sample_annotation.loc[
    sample_annotation["sample_type_code"].isin(tumor_codes),
    "sample_id"
].tolist()

tumor_expression = expression[tumor_samples].copy()

print("\nЧисло опухолевых образцов:", len(tumor_samples))

threshold = 1.0
minimum_fraction = 0.10
minimum_samples = int(len(tumor_samples) * minimum_fraction)

expressed_counts = (tumor_expression > threshold).sum(axis=1)

filtered_expression = tumor_expression[
    expressed_counts >= minimum_samples
].copy()

print("Порог экспрессии:", threshold)
print("Минимальное число образцов:", minimum_samples)
print("Псевдогенов до фильтрации:", tumor_expression.shape[0])
print("Псевдогенов после фильтрации:",
      filtered_expression.shape[0])

filtered_path = (
    processed_dir
    / "skcm_tumor_expressed_processed_pseudogenes.tsv.gz"
)

filtered_expression.to_csv(
    filtered_path,
    sep="\t",
    compression="gzip"
)

filtering_table = pd.DataFrame({
    "pseudogene": tumor_expression.index,
    "n_samples_expression_above_1": expressed_counts,
    "fraction_samples_expression_above_1":
        expressed_counts / len(tumor_samples),
    "mean_expression": tumor_expression.mean(axis=1),
    "median_expression": tumor_expression.median(axis=1)
}).sort_values(
    "mean_expression",
    ascending=False
)

filtering_table_path = (
    tables_dir
    / "pseudogene_expression_filtering.tsv"
)

filtering_table.to_csv(
    filtering_table_path,
    sep="\t",
    index=False
)

print("\nСохранены файлы:")
print(sample_annotation_path)
print(filtered_path)
print(filtering_table_path)

print("\nПсевдогены с наибольшей средней экспрессией:")
print(filtering_table.head(10).to_string(index=False))