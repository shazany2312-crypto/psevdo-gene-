from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

PROJECT_DIR = Path(__file__).resolve().parents[1]

input_path = (
    PROJECT_DIR
    / "results"
    / "tables"
    / "anxa2p1_mir186_potential_cerna_targets.tsv"
)

figures_dir = PROJECT_DIR / "figures"
figures_dir.mkdir(parents=True, exist_ok=True)

candidates = pd.read_csv(
    input_path,
    sep="\t"
)

G = nx.DiGraph()

pseudogene = "ANXA2P1"
mirna = "hsa-miR-186-5p"

G.add_node(pseudogene, node_type="pseudogene")
G.add_node(mirna, node_type="mirna")

for _, row in candidates.iterrows():
    target = row["target_gene"]

    G.add_node(target, node_type="target")

    G.add_edge(
        pseudogene,
        target,
        relation="positive",
        correlation=row["r_pseudogene_target"]
    )

    G.add_edge(
        mirna,
        target,
        relation="negative",
        correlation=row["r_mirna_target"]
    )

pos = {}

pos[pseudogene] = (-1.4, 0)
pos[mirna] = (0, 0)

targets = candidates["target_gene"].tolist()

for i, target in enumerate(targets):
    y = 4 - i * 0.6
    pos[target] = (2.1, y)

node_colors = []

for node in G.nodes:
    node_type = G.nodes[node]["node_type"]

    if node_type == "pseudogene":
        node_colors.append("#4C72B0")
    elif node_type == "mirna":
        node_colors.append("#C44E52")
    else:
        node_colors.append("#55A868")

plt.figure(figsize=(14, 10))

nx.draw_networkx_nodes(
    G,
    pos,
    node_color=node_colors,
    node_size=2200,
    edgecolors="black",
    linewidths=0.8
)

positive_edges = [
    (u, v)
    for u, v, data in G.edges(data=True)
    if data["relation"] == "positive"
]

negative_edges = [
    (u, v)
    for u, v, data in G.edges(data=True)
    if data["relation"] == "negative"
]

nx.draw_networkx_edges(
    G,
    pos,
    edgelist=positive_edges,
    edge_color="#4C72B0",
    width=1.8,
    arrows=True,
    arrowsize=18
)

nx.draw_networkx_edges(
    G,
    pos,
    edgelist=negative_edges,
    edge_color="#C44E52",
    style="dashed",
    width=1.8,
    arrows=True,
    arrowsize=18
)

nx.draw_networkx_labels(
    G,
    pos,
    font_size=9,
    font_weight="bold"
)

edge_labels = {}

for u, v, data in G.edges(data=True):
    edge_labels[(u, v)] = f"r = {data['correlation']:.2f}"

nx.draw_networkx_edge_labels(
    G,
    pos,
    edge_labels=edge_labels,
    font_size=8,
    label_pos=0.55
)

plt.title(
    "Потенциальная ceRNA-сеть ANXA2P1 — hsa-miR-186-5p — мРНК-мишени",
    fontsize=14
)

plt.text(
    -1.4,
    -1.2,
    "Синий: псевдоген\nКрасный: микроРНК\n"
    "Зелёный: мРНК-мишень\n"
    "Сплошная линия: положительная корреляция\n"
    "Пунктир: отрицательная корреляция",
    fontsize=10,
    bbox={
        "facecolor": "white",
        "edgecolor": "gray",
        "boxstyle": "round,pad=0.5"
    }
)

plt.axis("off")
plt.tight_layout()

output_path = (
    figures_dir
    / "potential_cerna_network_anxa2p1_mir186.png"
)

plt.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Число мРНК-кандидатов:", len(candidates))
print("Сеть сохранена:")
print(output_path)