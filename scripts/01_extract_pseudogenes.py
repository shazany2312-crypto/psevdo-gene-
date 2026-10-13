from pathlib import Path
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]

gtf_path = PROJECT_DIR / "data" / "gencode.v50.annotation.gtf.gz"
result_dir = PROJECT_DIR / "results" / "tables"
result_dir.mkdir(parents=True, exist_ok=True)

gtf_columns = [
    "chromosome",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute"
]

gtf = pd.read_csv(
    gtf_path,
    sep="\t",
    comment="#",
    header=None,
    names=gtf_columns
)

genes = gtf[gtf["feature"] == "gene"].copy()

genes["gene_id"] = genes["attribute"].str.extract(
    r'gene_id "([^"]+)"'
)

genes["gene_name"] = genes["attribute"].str.extract(
    r'gene_name "([^"]+)"'
)

genes["gene_type"] = genes["attribute"].str.extract(
    r'gene_type "([^"]+)"'
)

processed_pseudogenes = genes[
    genes["gene_type"] == "processed_pseudogene"
].copy()

columns_to_save = [
    "gene_id",
    "gene_name",
    "chromosome",
    "start",
    "end",
    "strand",
    "gene_type"
]

output_path = result_dir / "processed_pseudogenes.tsv"

processed_pseudogenes[columns_to_save].to_csv(
    output_path,
    sep="\t",
    index=False
)

print(f"Найдено псевдогенов: {len(processed_pseudogenes)}")
print(f"Результат сохранён: {output_path}")