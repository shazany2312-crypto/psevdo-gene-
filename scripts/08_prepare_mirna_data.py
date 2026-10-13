from pathlib import Path
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]

mirna_path = (
    PROJECT_DIR
    / "data"
    / "raw"
    / "TCGA_SKCM_sampleMap_miRNA_HiSeq_gene.gz"
)

pseudogene_path = (
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

processed_dir = PROJECT_DIR / "data" / "processed"
tables_dir = PROJECT_DIR / "results" / "tables"

print("1. Загружаем матрицу микроРНК...")

mirna = pd.read_csv(
    mirna_path,
    sep="\t",
    compression="gzip"
)

print("Размер матрицы микроРНК:", mirna.shape)
print("Первые 10 столбцов:")
print(mirna.columns[:10].tolist())
print("Первые 10 идентификаторов микроРНК:")
print(mirna.iloc[:10, 0].tolist())

mirna = mirna.rename(
    columns={mirna.columns[0]: "mirna_id"}
)

mirna["mirna_id"] = mirna["mirna_id"].astype(str).str.strip()

mirna = mirna.set_index("mirna_id")

print("\n2. Загружаем матрицу псевдогенов...")

pseudogene_expression = pd.read_csv(
    pseudogene_path,
    sep="\t",
    compression="gzip",
    index_col=0
)

annotation = pd.read_csv(
    annotation_path,
    sep="\t",
    dtype={"sample_type_code": str}
)

annotation["sample_type_code"] = (
    annotation["sample_type_code"]
    .str.zfill(2)
)

print("\nКоды типов образцов:")
print(annotation["sample_type_code"].value_counts())

tumor_samples = annotation.loc[
    annotation["sample_type_code"].isin(["01", "06", "07"]),
    "sample_id"
].tolist()

common_samples = sorted(
    set(mirna.columns)
    & set(pseudogene_expression.columns)
    & set(tumor_samples)
)

print("Опухолевых образцов в матрице псевдогенов:",
      len(tumor_samples))
print("Образцов в матрице микроРНК:",
      mirna.shape[1])
print("Общих опухолевых образцов:", len(common_samples))

mirna_common = mirna[common_samples].copy()
pseudogene_common = pseudogene_expression[common_samples].copy()

mirna_output_path = (
    processed_dir
    / "skcm_mirna_common_samples.tsv.gz"
)

pseudogene_output_path = (
    processed_dir
    / "skcm_pseudogene_common_samples.tsv.gz"
)

mirna_common.to_csv(
    mirna_output_path,
    sep="\t",
    compression="gzip"
)

pseudogene_common.to_csv(
    pseudogene_output_path,
    sep="\t",
    compression="gzip"
)

common_sample_table = annotation[
    annotation["sample_id"].isin(common_samples)
].copy()

common_sample_table_path = (
    tables_dir
    / "skcm_common_samples_mirna_pseudogenes.tsv"
)

common_sample_table.to_csv(
    common_sample_table_path,
    sep="\t",
    index=False
)

print("\nСохранены файлы:")
print(mirna_output_path)
print(pseudogene_output_path)
print(common_sample_table_path)