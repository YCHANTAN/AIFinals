import pandas as pd
import matplotlib.pyplot as plt
import os

# -----------------------
# Load results
# -----------------------
df = pd.read_csv("results/results.csv")

# Drop invalid rows
df = df[df["tanimoto"].notna()]

# -----------------------
# Plot histogram
# -----------------------
plt.figure(figsize=(6, 4))
plt.hist(df["tanimoto"], bins=20)
plt.xlabel("Tanimoto similarity")
plt.ylabel("Count")
plt.title("Reconstruction Similarity Distribution")

# -----------------------
# Save figure
# -----------------------
os.makedirs("figures", exist_ok=True)
plt.savefig("figures/tanimoto_histogram.png", dpi=300, bbox_inches="tight")
plt.close()

print("Saved figures/tanimoto_histogram.png")
