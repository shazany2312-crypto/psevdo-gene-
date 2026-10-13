from pathlib import Path
import pandas as pd

DATA_DIR = Path(
    r"C:\Users\Пользователь\Desktop\melanoma_pseudogenes\melanoma_pseudogenes\data\raw"
)

print("Папка с данными:")
print(DATA_DIR)

print("\nВсе файлы в этой папке:")
for file in DATA_DIR.iterdir():
    print(file.name)

mirna_files = list(DATA_DIR.glob("*miRNA*"))

print("\nНайденные miRNA-файлы:")
for file in mirna_files:
    print(file.name)

if len(mirna_files) == 0:
    raise FileNotFoundError(
        "Файл с miRNA не найден. Проверь, что он действительно находится в указанной папке."
    )

mirna_path = mirna_files[0]

mirna = pd.read_csv(
    mirna_path,
    sep="\t",
    compression="gzip",
    low_memory=False
)

print("\nИспользуемый файл:")
print(mirna_path.name)

print("\nРазмер таблицы:")
print(mirna.shape)

print("\nПервые 10 строк и первые 5 столбцов:")
print(mirna.iloc[:10, :5])

print("\nПервые 10 названий столбцов:")
print(mirna.columns[:10].tolist())