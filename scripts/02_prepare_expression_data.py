from pathlib import Path
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]

raw_dir = PROJECT_DIR / "data" / "raw"

phenotype_path = raw_dir / "TcgaTargetGTEX_phenotype.txt.gz"

phenotype = pd.read_csv(
    phenotype_path,
    sep="\t",
    low_memory=False
)

print("Размер таблицы метаданных:", phenotype.shape)
print("\nНазвания столбцов:")
print(phenotype.columns.tolist())

print("\nПервые 5 строк:")
print(phenotype.head())