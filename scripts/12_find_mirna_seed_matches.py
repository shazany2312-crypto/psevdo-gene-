from pathlib import Path
import pandas as pd


def reverse_complement(sequence):
    complement = str.maketrans(
        "ACGTUacgtu",
        "TGCAAtgcaa"
    )

    return sequence.translate(complement)[::-1].upper()


def read_fasta(fasta_path):
    sequences = {}
    current_name = None
    current_sequence = []

    with open(fasta_path, "r", encoding="utf-8") as fasta_file:
        for line in fasta_file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if current_name is not None:
                    sequences[current_name] = "".join(
                        current_sequence
                    ).upper()

                current_name = line[1:].split()[0]
                current_sequence = []
            else:
                current_sequence.append(line)

        if current_name is not None:
            sequences[current_name] = "".join(
                current_sequence
            ).upper()

    return sequences


PROJECT_DIR = Path(__file__).resolve().parents[1]

candidate_path = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "anxa2_pseudogene_negative_mirna_candidates_named.tsv"
)

mirna_fasta_path = (
    PROJECT_DIR
    / "data"
    / "raw"
    / "mature.fa"
)

pseudogene_fasta_path = (
    PROJECT_DIR
    / "data"
    / "raw"
    / "anxa2p_pseudogenes_GRCh38.fasta"
)

tables_dir = PROJECT_DIR / "results" / "tables"

candidates = pd.read_csv(
    candidate_path,
    sep="\t"
)

mirna_sequences = read_fasta(mirna_fasta_path)
pseudogene_sequences = read_fasta(pseudogene_fasta_path)

human_mirna_sequences = {}

with open(mirna_fasta_path, "r", encoding="utf-8") as fasta_file:
    current_header = None
    current_sequence = []

    for line in fasta_file:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            if current_header is not None:
                parts = current_header[1:].split()
                mirna_name = parts[0]
                mirna_id = next(
                    (
                        part
                        for part in parts
                        if part.startswith("MIMAT")
                    ),
                    None
                )

                if (
                    mirna_name.startswith("hsa-")
                    and mirna_id is not None
                ):
                    human_mirna_sequences[mirna_id] = {
                        "mirna_name": mirna_name,
                        "sequence": "".join(current_sequence)
                        .upper()
                        .replace("U", "T")
                    }

            current_header = line
            current_sequence = []
        else:
            current_sequence.append(line)

    if current_header is not None:
        parts = current_header[1:].split()
        mirna_name = parts[0]
        mirna_id = next(
            (
                part
                for part in parts
                if part.startswith("MIMAT")
            ),
            None
        )

        if mirna_name.startswith("hsa-") and mirna_id is not None:
            human_mirna_sequences[mirna_id] = {
                "mirna_name": mirna_name,
                "sequence": "".join(current_sequence)
                .upper()
                .replace("U", "T")
            }

print("Корреляционных кандидатов:", len(candidates))
print("Последовательностей псевдогенов:", len(pseudogene_sequences))

matches = []

for _, row in candidates.iterrows():
    pseudogene = row["pseudogene"]
    mirna_id = row["mirna_id"]

    if mirna_id not in human_mirna_sequences:
        continue

    mirna_name = human_mirna_sequences[mirna_id]["mirna_name"]
    mirna_sequence = human_mirna_sequences[mirna_id]["sequence"]

    if len(mirna_sequence) < 8:
        continue

    seed = mirna_sequence[1:8]
    target_seed = reverse_complement(seed)

    pseudogene_sequence = pseudogene_sequences[pseudogene]

    start_position = 0

    while True:
        position = pseudogene_sequence.find(
            target_seed,
            start_position
        )

        if position == -1:
            break

        matches.append({
            "pseudogene": pseudogene,
            "mirna_id": mirna_id,
            "mirna_name": mirna_name,
            "mirna_sequence_5_to_3": mirna_sequence,
            "seed_2_to_8": seed,
            "target_seed_match_5_to_3": target_seed,
            "match_start_1_based": position + 1,
            "match_end_1_based": position + len(target_seed),
            "spearman_r": row["spearman_r"],
            "FDR": row["FDR"]
        })

        start_position = position + 1

matches = pd.DataFrame(matches)

output_path = (
    tables_dir
    / "anxa2p_mirna_seed_match_candidates.tsv"
)

matches.to_csv(
    output_path,
    sep="\t",
    index=False
)

print("\nНайдено точных 7-нуклеотидных seed-match:",
      len(matches))

if len(matches) > 0:
    print("\nНайденные кандидаты:")
    print(
        matches[
            [
                "pseudogene",
                "mirna_name",
                "seed_2_to_8",
                "target_seed_match_5_to_3",
                "match_start_1_based",
                "match_end_1_based",
                "spearman_r",
                "FDR"
            ]
        ].to_string(index=False)
    )
else:
    print(
        "\nСреди 24 корреляционных кандидатов "
        "точных 7-нуклеотидных seed-match не найдено."
    )

print("\nТаблица сохранена:")
print(output_path)