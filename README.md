# BiCLUM_test

This repository contains the code and resources for reproducing the results and evaluations presented in the paper **BiCLUM: Bilateral Contrastive Learning for Unpaired Single-Cell Multi-Omics Integration**. It includes implementations of the comparison methods, quantitative evaluation metrics, and visualization scripts.

---

## Repository Structure

The repository is organized into three main folders:

1. **`compared_methods`**  
   This folder contains the implementation of all comparison methods used in the paper. Each method includes a script for running the analysis on the PBMC dataset (e.g., `PBMC.py`).

2. **`eva`**  
   This folder includes the scripts for evaluating the results using the quantitative metrics described in the paper. For example, `eva_PBMC.py` computes all metrics for the PBMC dataset.

3. **`vis`**  
   This folder contains scripts for generating visualizations of the results. For instance, `vis_pbmc.py` produces all figures related to the PBMC dataset included in the paper.


---

## Included Comparison Methods  

Below is a list of the comparison methods evaluated in this study, along with their respective GitHub repositories:  

- **[GLUE](https://github.com/gao-lab/GLUE)**  
- **[JointMDS](https://github.com/BorgwardtLab/JointMDS)**  
- **[LIGER](https://github.com/welch-lab/liger)**  
- **[MMD-MA](https://bitbucket.org/noblelab/2019_mmd_wabi/src/master/)**  
- **[MultiMAP](https://github.com/Teichlab/MultiMAP)**  
- **[Pamona](https://github.com/caokai1073/Pamona)**  
- **[SCOT](https://github.com/rsinghlab/SCOT)**  
- **[Seurat](https://satijalab.org/seurat/)**  
- **[UnionCom](https://github.com/caokai1073/UnionCom)**  
- **[bindSC](https://github.com/KChen-lab/bindSC)**  
- **[scDART](https://github.com/PeterZZQ/scDART)**  
- **[scTopoGAN](https://github.com/AkashCiel/scTopoGAN)**  
- **[uniPort](https://github.com/caokai1073/uniPort)**  

---

## How to Reproduce Results

### Step 1: Run Comparison Methods
1. Navigate to the `compared_methods` folder.
2. Run the `PBMC.py` script for each comparison method.  
   Example:
   ```bash
   python compared_methods/method_name/PBMC.py
   ```
   - Replace `method_name` with the desired method (e.g., `Seurat`, `GLUE`).
   - The results for each method will be saved in a `results` folder specific to that method.

### Step 2: Evaluate Results
1. Navigate to the `eva` folder.
2. Run the evaluation script for the PBMC dataset:  
   ```bash
   python eva/eva_PBMC.py
   ```
   - This script will compute the quantitative metrics, such as omics mixing, biology conservation, transfer accuracy, and FOSCTTM scores.

### Step 3: Generate Visualizations
1. Navigate to the `vis` folder.
2. Run the visualization script for the PBMC dataset:  
   ```bash
   python vis/vis_pbmc.py
   ```
   - This script will generate all visualizations, including UMAP plots, PAGA graphs, and comparison figures, and save them in a designated `plots` folder.

---

## Notes
- **File Outputs**: Each step will generate outputs in corresponding subdirectories:
  - Comparison results: `results/method_name/`
  - Evaluation metrics: `eva/`
  - Visualizations: `vis/`
- **PBMC Example**: While the instructions focus on PBMC data, similar steps can be followed for other datasets by replacing `PBMC` with the corresponding dataset name in the scripts.

---  

For further details, refer to the documentation and supplementary materials in the repository.
