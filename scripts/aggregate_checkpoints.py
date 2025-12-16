import pandas as pd
import glob

rows = []

for file in glob.glob("results/checkpoints/results_*.csv"):
    step = file.split("_")[-1].replace(".csv", "")
    df = pd.read_csv(file)

    rows.append({
        "step": int(step),
        "exact_match_rate": df["exact_match"].mean(),
        "mean_tanimoto": df["tanimoto"].mean()
    })

summary = pd.DataFrame(rows).sort_values("step")
summary.to_csv("results/checkpoint_summary.csv", index=False)
print(summary)
