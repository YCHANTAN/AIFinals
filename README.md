# Measuring Information Loss in Hierarchical Graph-to-Graph Molecular Autoencoders

## 📖 Exam Context & Objective
This repository contains the solution for the **CS 3rd Year Take-Home Exam**.
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
