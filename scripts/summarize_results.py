import pandas as pd
import os

df = pd.read_csv("results/results.csv")

summary = {
    "total_samples": len(df),
    "valid_rate": df["valid_out"].mean(),
    "exact_match_rate": df["exact_match"].mean(),
    "mean_tanimoto": df["tanimoto"].mean(),
    "mean_delta_atoms": df["delta_atoms"].mean(),
    "mean_delta_bonds": df["delta_bonds"].mean()
}

summary_df = pd.DataFrame([summary])

os.makedirs("results", exist_ok=True)
summary_df.to_csv("results/summary_metrics.csv", index=False)

print("Saved results/summary_metrics.csv")
print(summary_df)
