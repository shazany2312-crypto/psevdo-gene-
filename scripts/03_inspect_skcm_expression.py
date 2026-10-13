from pathlib import Path
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]

expression_path = (
    PROJECT_DIR
    / "data"
    / "raw"
    / "TCGA_SKCM_sampleMap_HiSeqV2.gz"
)

expression = pd.read_csv(
    expression_path,
    sep="\t",
    compression="gzip",
    nrows=10
)

print("Размер прочитанного фрагмента:", expression.shape)

print("\nНазвания первых 10 столбцов:")
print(expression.columns[:10].tolist())

print("\nПервые 10 строк:")
print(expression.head(10).to_string())