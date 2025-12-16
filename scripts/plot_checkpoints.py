import pandas as pd
import matplotlib.pyplot as plt
import os

df = pd.read_csv("results/checkpoint_summary.csv")

os.makedirs("figures", exist_ok=True)

plt.figure()
plt.plot(df["step"], df["exact_match_rate"], marker="o")
plt.xlabel("Training step")
plt.ylabel("Exact match accuracy")
plt.title("Exact Match vs Training Step")
plt.savefig("figures/exact_match_vs_step.png", dpi=300)
plt.close()

plt.figure()
plt.plot(df["step"], df["mean_tanimoto"], marker="o")
plt.xlabel("Training step")
plt.ylabel("Mean Tanimoto")
plt.title("Tanimoto Similarity vs Training Step")
plt.savefig("figures/tanimoto_vs_step.png", dpi=300)
plt.close()
