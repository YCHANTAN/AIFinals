from chembl_webresource_client.new_client import new_client
import pandas as pd
import os
import warnings

warnings.filterwarnings("ignore", category=UserWarning)

# -----------------------
# STEP 1: Choose target
# -----------------------
TARGET_ID = "CHEMBL203"  # Canonical human EGFR
print("Using target:", TARGET_ID)

# -----------------------
# STEP 2: Query activities
# -----------------------
activity = new_client.activity

results = activity.filter(
    target_chembl_id=TARGET_ID,
    standard_type="IC50"
).only(
    "molecule_chembl_id",
    "canonical_smiles",
    "standard_value",
    "standard_units"
)

# -----------------------
# STEP 3: Collect data (LIMITED)
# -----------------------
data = []
MAX_RECORDS = 1000  # SAFE LIMIT

for i, r in enumerate(results):
    if i >= MAX_RECORDS:
        break

    if r.get("canonical_smiles") and r.get("standard_value"):
        data.append({
            "molecule_id": r["molecule_chembl_id"],
            "smiles": r["canonical_smiles"],
            "ic50": r["standard_value"],
            "units": r["standard_units"]
        })

# -----------------------
# STEP 4: Save CSV
# -----------------------
df = pd.DataFrame(data)

os.makedirs("data", exist_ok=True)
df.to_csv("data/chembl_egfr_ic50.csv", index=False)

print("Saved data/chembl_egfr_ic50.csv")
print("Total molecules:", len(df))
