import sys
if sys.platform == "win32":
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", errors="replace")

import argparse, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog("rdApp.*")

FAMILY_COLORS = {"benz": "#4C72B0", "naph": "#DD8452", "ind": "#55A868",
                 "quin": "#C44E52", "pyr": "#8172B2", "bzim": "#937860", "other": "#808080"}

# Ionizable group SMARTS + literature pKa defaults
IONIZABLE_GROUPS = {
    "aliphatic_amine": ("[C][NH2,NH1;!$(NC=O);!$(Nc):1]", 10.0, "base"),
    "aromatic_amine":  ("[c][NH2,NH1;!$(NC=O):1]",         4.5,  "base"),
    "carboxylic_acid": ("[C:1](=O)[OH1]",                  4.5,  "acid"),
    "phenol":          ("[c:1][OH1]",                       10.0, "acid"),
    "sulfonamide_nh":  ("[S](=O)(=O)[NH1:1]",              10.5, "acid"),
    "imidazole":       ("[nH:1]1ccnc1",                    6.0,  "base"),
}

def load_compounds(path):
    df = pd.read_csv(path)
    records, n_bad = [], 0
    normalizer = rdMolStandardize.Normalizer()
    largest_frag = rdMolStandardize.LargestFragmentChooser()
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["smiles"]))
        if mol is None:
            n_bad += 1
            continue
        mol = normalizer.normalize(mol)
        mol = largest_frag.choose(mol)
        fam = str(row["compound_name"]).split("_")[0]
        records.append({
            "compound_name": str(row["compound_name"]),
            "family": fam if fam in FAMILY_COLORS else "other",
            "pic50": float(row.get("pic50", float("nan"))),
            "mol": mol,
        })
    print(f"  {len(records)} valid ({n_bad} skipped)")
    return pd.DataFrame(records)

def compute_charge_at_ph(mol, ph=7.4):
    """Estimate formal charge at given pH using pKa rules."""
    patterns = {name: Chem.MolFromSmarts(sma) for name, (sma, pka, kind) in IONIZABLE_GROUPS.items()}
    total_charge = 0
    n_ionizable = 0
    group_detail = []

    for name, (sma, pka, kind) in IONIZABLE_GROUPS.items():
        pat = patterns[name]
        if pat is None:
            continue
        matches = mol.GetSubstructMatches(pat)
        for _ in matches:
            n_ionizable += 1
            if kind == "base":
                # base: pKa of conjugate acid; protonated (charge +1) when pH < pKa
                charged = ph < pka
                total_charge += 1 if charged else 0
            else:
                # acid: deprotonated (charge -1) when pH > pKa
                charged = ph > pka
                total_charge += -1 if charged else 0
            group_detail.append(name)

    return total_charge, n_ionizable

def apply_protonation(df, ph=7.4):
    rows = []
    for _, row in df.iterrows():
        charge, n_ionizable = compute_charge_at_ph(row["mol"], ph)
        rows.append({
            "compound_name": row["compound_name"],
            "family": row["family"],
            "pic50": row["pic50"],
            "n_ionizable_groups": n_ionizable,
            "formal_charge_ph74": charge,
            "dominant_form": f"{'+' if charge > 0 else ''}{charge}",
        })
    return pd.DataFrame(rows)

def plot_charge_distribution(result_df, output_path):
    families = [f for f in FAMILY_COLORS if f in result_df["family"].values]
    charge_values = sorted(result_df["formal_charge_ph74"].unique())
    charge_labels = [f"{'+' if c > 0 else ''}{c}" for c in charge_values]

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(families))
    bar_width = 0.8 / max(len(charge_values), 1)
    bottom = np.zeros(len(families))

    charge_colors = {-2: "#d73027", -1: "#fc8d59", 0: "#ffffbf", 1: "#91bfdb", 2: "#4575b4"}

    counts_matrix = {}
    for c in charge_values:
        counts_matrix[c] = []
        for fam in families:
            sub = result_df[result_df["family"] == fam]
            counts_matrix[c].append((sub["formal_charge_ph74"] == c).sum())

    for c in charge_values:
        counts = np.array(counts_matrix[c])
        color = charge_colors.get(c, "#808080")
        label = f"{'+' if c > 0 else ''}{c}"
        ax.bar(x, counts, bottom=bottom, color=color, label=label, edgecolor="white", linewidth=0.5)
        bottom += counts

    ax.set_xticks(x)
    ax.set_xticklabels(families, fontsize=10)
    ax.set_xlabel("Scaffold Family", fontsize=11)
    ax.set_ylabel("Number of Compounds", fontsize=11)
    ax.set_title(f"Charge State Distribution at pH 7.4", fontsize=13, fontweight="bold")
    ax.legend(title="Formal Charge", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--ph", type=float, default=7.4)
    parser.add_argument("--output-dir", default="output")
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"\nLoading: {args.input}")
    df = load_compounds(args.input)
    result_df = apply_protonation(df, ph=args.ph)

    csv_path = os.path.join(args.output_dir, "protonation_states.csv")
    result_df.to_csv(csv_path, index=False)
    print(f"Saved: {csv_path}")

    plot_charge_distribution(result_df, os.path.join(args.output_dir, "charge_distribution.png"))
    print(f"Saved: {args.output_dir}/charge_distribution.png")

    print(f"\n--- Charge distribution at pH {args.ph} ---")
    charge_summary = result_df.groupby("formal_charge_ph74")["compound_name"].count().rename("count")
    print(charge_summary.to_string())

    print("\n--- Ionizable groups summary by family ---")
    summary = result_df.groupby("family")["n_ionizable_groups"].agg(["mean", "min", "max"]).round(2)
    print(summary.to_string())
    print("\nDone.")

if __name__ == "__main__":
    main()
