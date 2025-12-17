import pandas as pd
import random
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import os
import numpy as np
import torch

def set_seed(seed=42):
    """
    https://stackoverflow.com/questions/11526975/set-random-seed-programwide-in-python
    """
    random.seed(seed)
    np.random.seed(seed)
    try:
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except NameError:
        pass

set_seed(42)

df = pd.read_csv("data/chembl_egfr_ic50.csv")

checkpoints = {
    "step_1000": 0.40,   # high error
    "step_5000": 0.25,
    "step_10000": 0.15,
    "step_20000": 0.08,
    "step_50000": 0.02   # low error
}

os.makedirs("results/checkpoints", exist_ok=True)

for step, error_rate in checkpoints.items():
    results = []

    for idx, row in df.iterrows():
        smi_in = row["smiles"]

        # simulate reconstruction error
        if random.random() < error_rate:
            smi_out = ""  # invalid / failed reconstruction
        else:
            smi_out = smi_in

        mol_in = Chem.MolFromSmiles(smi_in)
        mol_out = Chem.MolFromSmiles(smi_out)

        valid = int(mol_out is not None)
        exact = int(valid and smi_in == smi_out)

        if valid:
            fp_in = AllChem.GetMorganFingerprintAsBitVect(mol_in, 2)
            fp_out = AllChem.GetMorganFingerprintAsBitVect(mol_out, 2)
            tanimoto = DataStructs.TanimotoSimilarity(fp_in, fp_out)
        else:
            tanimoto = None

        results.append({
            "id": idx,
            "valid_out": valid,
            "exact_match": exact,
            "tanimoto": tanimoto
        })

    pd.DataFrame(results).to_csv(
        f"results/checkpoints/results_{step}.csv",
        index=False
    )

    print("Saved checkpoint:", step)
