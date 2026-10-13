from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

DATA_DIR = Path(
    r"C:\Users\Пользователь\Desktop\melanoma_pseudogenes\melanoma_pseudogenes\data\raw"
)

RESULTS_DIR = Path(
    r"C:\Users\Пользователь\Desktop\melanoma_pseudogenes\melanoma_pseudogenes\results"
)

FIGURES_DIR = Path(
    r"C:\Users\Пользователь\Desktop\melanoma_pseudogenes\melanoma_pseudogenes\figures"
)

RESULTS_DIR.mkdir(exist_ok=True)
FIGURES_DIR.mkdir(exist_ok=True)

mrna_path = DATA_DIR / "TCGA_SKCM_sampleMap_HiSeqV2.gz"
mirna_path = DATA_DIR / "TCGA_SKCM_sampleMap_miRNA_HiSeq_gene.gz"

PSEUDOGENE = "ANXA2P1"
MIRNA_ID = "MIMAT0000456"
MIRNA_NAME = "hsa-miR-186-5p"

mrna = pd.read_csv(
    mrna_path,
    sep="\t",
    compression="gzip",
    index_col=0
)

mirna = pd.read_csv(
    mirna_path,
    sep="\t",
    compression="gzip",
    index_col=0
)

print("Размер матрицы генов:", mrna.shape)
print("Размер матрицы miRNA:", mirna.shape)

print("\nПроверяем наличие ANXA2P1:")
print(PSEUDOGENE in mrna.index)

print("\nПроверяем наличие MIMAT0000456:")
print(MIRNA_ID in mirna.index)

if PSEUDOGENE not in mrna.index:
    matches = [
        index for index in mrna.index
        if PSEUDOGENE.lower() in str(index).lower()
    ]
    print("\nПохожие строки в матрице генов:")
    print(matches[:20])

    raise KeyError(
        "ANXA2P1 не найден. Посмотри похожие названия выше."
    )

if MIRNA_ID not in mirna.index:
    raise KeyError(
        "MIMAT0000456 не найден в miRNA-матрице."
    )

anxa2p1 = mrna.loc[PSEUDOGENE]
mir186 = mirna.loc[MIRNA_ID]

common_samples = anxa2p1.index.intersection(mir186.index)

anxa2p1 = anxa2p1[common_samples]
mir186 = mir186[common_samples]

data = pd.DataFrame({
    "ANXA2P1": anxa2p1,
    MIRNA_NAME: mir186
}).dropna()

rho, p_value = spearmanr(
    data["ANXA2P1"],
    data[MIRNA_NAME]
)

result = pd.DataFrame([{
    "pseudogene": PSEUDOGENE,
    "miRNA": MIRNA_NAME,
    "miRNA_id": MIRNA_ID,
    "n_common_samples_before_na_removal": len(common_samples),
    "n_samples_used": len(data),
    "spearman_rho": rho,
    "p_value": p_value
}])

result_path = RESULTS_DIR / "anxa2p1_mir186_spearman.tsv"

result.to_csv(
    result_path,
    sep="\t",
    index=False
)

plt.figure(figsize=(7, 6))

sns.regplot(
    data=data,
    x="ANXA2P1",
    y=MIRNA_NAME,
    scatter_kws={"alpha": 0.65, "s": 30},
    line_kws={"color": "crimson"}
)

plt.title(
    f"ANXA2P1 и {MIRNA_NAME}\n"
    f"Spearman ρ = {rho:.3f}, p = {p_value:.3e}, n = {len(data)}"
)
plt.xlabel("Экспрессия ANXA2P1")
plt.ylabel(f"Экспрессия {MIRNA_NAME}")
plt.tight_layout()

figure_path = FIGURES_DIR / "anxa2p1_mir186_correlation.png"

plt.savefig(figure_path, dpi=300)
plt.show()

print("\nЧисло общих образцов:", len(common_samples))
print("Число образцов после удаления пропусков:", len(data))
print(f"Корреляция Спирмена: rho = {rho:.4f}")
print(f"p-value: {p_value:.6g}")

print("\nТаблица сохранена:")
print(result_path)

print("\nГрафик сохранён:")
print(figure_path)