from pathlib import Path

DATA_DIR = Path(
    r"C:\Users\Пользователь\Desktop\melanoma_pseudogenes\melanoma_pseudogenes\data\raw"
)

mature_fasta = DATA_DIR / "mature.fa"

target_name = "hsa-miR-186-5p"

with open(mature_fasta, "r", encoding="utf-8") as file:
    for line in file:
        if line.startswith(">") and target_name in line:
            print("Найдена строка:")
            print(line.strip())