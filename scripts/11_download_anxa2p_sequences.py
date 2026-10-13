from pathlib import Path
import time
import pandas as pd
import requests

PROJECT_DIR = Path(__file__).resolve().parents[1]

annotation_path = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "processed_pseudogenes.tsv"
)

raw_dir = PROJECT_DIR / "data" / "raw"
tables_dir = PROJECT_DIR / "results" / "tables"

pseudogenes = pd.read_csv(
    annotation_path,
    sep="\t"
)

candidates = ["ANXA2P1", "ANXA2P2", "ANXA2P3"]

selected = pseudogenes[
    pseudogenes["gene_name"].isin(candidates)
].copy()

selected = selected[
    [
        "gene_id",
        "gene_name",
        "chromosome",
        "start",
        "end",
        "strand",
        "gene_type"
    ]
].sort_values("gene_name")

print("Найдены в аннотации:")
print(selected.to_string(index=False))

if len(selected) != 3:
    raise ValueError(
        "Не удалось найти все три псевдогена в таблице аннотации."
    )

headers = {
    "Content-Type": "text/plain",
    "User-Agent": "melanoma-pseudogene-coursework"
}

fasta_records = []
download_info = []

for _, row in selected.iterrows():
    chromosome = str(row["chromosome"]).replace("chr", "")

    strand_number = 1 if row["strand"] == "+" else -1

    region = (
        f"{chromosome}:"
        f"{row['start']}..{row['end']}:"
        f"{strand_number}"
    )

    url = (
        "https://rest.ensembl.org/"
        f"sequence/region/human/{region}"
    )

    print(f"\nСкачиваем {row['gene_name']}: {region}")

    response = requests.get(
        url,
        headers=headers,
        timeout=30
    )

    if response.status_code != 200:
        raise RuntimeError(
            f"Ошибка Ensembl для {row['gene_name']}: "
            f"{response.status_code}\n{response.text}"
        )

    sequence = (
        response.text
        .replace("\n", "")
        .strip()
        .upper()
    )

    if not sequence:
        raise ValueError(
            f"Ensembl вернул пустую последовательность "
            f"для {row['gene_name']}."
        )

    fasta_records.append(
        f">{row['gene_name']} "
        f"| {row['gene_id']} "
        f"| GRCh38 "
        f"| {region}\n"
        f"{sequence}\n"
    )

    download_info.append({
        "gene_id": row["gene_id"],
        "gene_name": row["gene_name"],
        "chromosome": row["chromosome"],
        "start": row["start"],
        "end": row["end"],
        "strand": row["strand"],
        "ensembl_region": region,
        "sequence_length": len(sequence)
    })

    print(f"Длина последовательности: {len(sequence)} нт")

    time.sleep(0.2)

fasta_path = raw_dir / "anxa2p_pseudogenes_GRCh38.fasta"

with open(fasta_path, "w", encoding="utf-8") as fasta_file:
    fasta_file.writelines(fasta_records)

info_path = tables_dir / "anxa2p_sequence_download_info.tsv"

pd.DataFrame(download_info).to_csv(
    info_path,
    sep="\t",
    index=False
)

print("\nГотово.")
print("FASTA-файл:", fasta_path)
print("Таблица координат и длин:", info_path)