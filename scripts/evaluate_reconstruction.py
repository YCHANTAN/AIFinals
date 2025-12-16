import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import os
gen = GetMorganGenerator(radius=2, fpSize=2048)

# -----------------------
# Load input molecules
# -----------------------
INPUT_CSV = "data/chembl_egfr_ic50.csv"
df = pd.read_csv(INPUT_CSV)

results = []

for idx, row in df.iterrows():
    smi_in = row["smiles"]

    # === SANITY RECONSTRUCTION ===
    smi_out = smi_in  # replace later with model output

    mol_in = Chem.MolFromSmiles(smi_in)
    mol_out = Chem.MolFromSmiles(smi_out)

    valid_out = int(mol_out is not None)

    if valid_out:
        exact = int(
            Chem.MolToSmiles(mol_in, canonical=True)
            == Chem.MolToSmiles(mol_out, canonical=True)
        )

        fp_in = gen.GetFingerprint(mol_in)
        fp_out = gen.GetFingerprint(mol_out)
        tanimoto = DataStructs.TanimotoSimilarity(fp_in, fp_out)

        delta_atoms = abs(mol_in.GetNumAtoms() - mol_out.GetNumAtoms())
        delta_bonds = abs(mol_in.GetNumBonds() - mol_out.GetNumBonds())
    else:
        exact = 0
        tanimoto = None
        delta_atoms = None
        delta_bonds = None

    results.append({
        "id": idx,
        "smiles_in": smi_in,
        "smiles_out": smi_out if valid_out else "",
        "valid_out": valid_out,
        "exact_match": exact,
        "tanimoto": tanimoto,
        "delta_atoms": delta_atoms,
        "delta_bonds": delta_bonds
    })

# -----------------------
# Save results
# -----------------------
os.makedirs("results", exist_ok=True)
out_path = "results/results.csv"
pd.DataFrame(results).to_csv(out_path, index=False)

print("Saved", out_path)
