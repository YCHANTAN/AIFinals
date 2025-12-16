import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import os

# -----------------------
# Load results
# -----------------------
df = pd.read_csv("results/results.csv")

# Keep only valid reconstructions
df = df[df["valid_out"] == 1]

# Take 10 examples
examples = df.head(10)

mols = []
legends = []

for _, row in examples.iterrows():
    mol_in = Chem.MolFromSmiles(row["smiles_in"])
    mol_out = Chem.MolFromSmiles(row["smiles_out"])

    mols.append(mol_in)
    legends.append("Input")

    mols.append(mol_out)
    legends.append("Reconstruction")

# -----------------------
# Draw molecules
# -----------------------
img = Draw.MolsToGridImage(
    mols,
    molsPerRow=2,
    subImgSize=(250, 250),
    legends=legends
)

# -----------------------
# Save image
# -----------------------
os.makedirs("figures", exist_ok=True)
out_path = "figures/qualitative_examples.png"
img.save(out_path)

print("Saved", out_path)
