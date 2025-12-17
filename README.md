# Measuring Information Loss in Hierarchical Graph-to-Graph Molecular Autoencoders

## 📖 Exam Context & Objective
This repository contains the solution for the **Information Loss in
Hierarchical Graph-to-Graph Molecular Autoencoders**.
The objective is to quantify information loss in the `hgraph2graph` Variational Autoencoder (VAE) by measuring reconstruction accuracy, structural similarity, and decoder dynamics.

**Core Tasks:**
1.  **Part A:** Quantify reconstruction loss (Exact Match, Tanimoto, Validity).
2.  **Part B:** Analyze checkpoint dynamics during training.
3.  **Part C:** Interpret decoder behavior (Hidden-State Probing or Partial Decoding).
4.  **Part D:** Explain the training mechanics.

---

## 🛠️ Environment Setup (Reproducibility)
To ensure the code runs end-to-end on a clean machine, please follow these steps.

**1. Prerequisites**
* Python 3.7+
* Anaconda or Miniconda

**2. Installation**
Create the environment and install dependencies (PyTorch, RDKit, NumPy, Pandas, Matplotlib).

**2. Clone Repo**
```bash
git clone https://github.com/wengong-jin/hgraph2graph.git


```bash
# Create the environment
conda create -n hgraph_exam python=3.8 -y
conda activate hgraph_exam

# Install RDKit (Required for chemical validity checks)
conda install -c conda-forge rdkit

# Install PyTorch (Adjust for your CUDA version)
conda install pytorch torchvision -c pytorch

# Install standard analysis tools
pip install pandas matplotlib numpy tqdm

## Reproducibility & Nondeterminism

To ensure the experiments are reproducible, all random seeds have been fixed to `42` using the `set_seed()` function provided in the scripts. This initializes the random number generators for:
* Python (`random`)
* NumPy (`np.random`)
* PyTorch (`torch.manual_seed`)

### Note on Nondeterminism
While all seeds are fixed, minor nondeterminism may still persist in the following areas:
1.  **GPU Floating Point Operations:** If training is run on a GPU, parallel computations in CUDA can introduce slight variations in loss values compared to CPU execution.
2.  **Hardware Differences:** Running this code on different architectures (e.g., Google Colab T4 vs. local RTX 3060) may result in negligible differences in the final trained weights.

The results submitted in `results.csv` were generated using [INSERT YOUR DEVICE HERE, e.g., Google Colab T4 GPU].
