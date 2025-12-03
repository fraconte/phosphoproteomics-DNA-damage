# phosphoproteomics-DNA-damage
# DNA Damage Response Phosphoproteomics

This repository contains analysis scripts and supporting code for a study investigating the dynamics of protein phosphorylation during the **DNA Damage Response (DDR)** across multiple genotoxic treatments.

## Background

The DNA damage response (DDR) is an intricate network of protein–protein interactions and signaling pathways activated by DNA lesions and genomic instability. A central regulatory mechanism of this response is **protein phosphorylation**, which coordinates DNA repair, activation of cell cycle checkpoints, and chromatin organization.

Despite its importance, the global behavior of the human phosphoproteome following exposure to diverse DNA damage–inducing agents remains insufficiently characterized. Understanding how phosphorylation sites are regulated in response to specific types of DNA damage is essential to define shared and pathway-specific DDR signaling signatures.

## Study Overview

In this project, we systematically profiled the cellular phosphoproteome after treatment with **eleven different DNA damage–inducing agents** encountered physiologically or used in cancer therapy. We investigated treatment-specific signatures, common responses, and the enrichment of regulated proximal phosphorylation sites in (phospho-clustering) in intrinsically disordered regions (IDRs).

## Repository Contents

This repository contains scripts used for:

- Phosphoproteomics data preprocessing and filtering
- Differential phosphorylation analysis across treatments
- Functional annotation and enrichment analyses
- Phosphosite clustering
- Data visualization and figure generation


## Requirements

Analyses were performed using **R** (≥4.2).  
Key R packages include:

- `clusterProfiler`
- `DOSE`
- `dplyr`
- `forcats`
- `ggplot2`
- `pheatmap`
- `stringr`
- `tidyr`

Additional dependencies are documented within individual scripts.

## Citation

If you use any part of this repository or data in your work, please cite the corresponding study:


## Contact

For questions or issues related to this repository, please contact:

**Francesca Conte**  
f.conte@imb-mainz.de

---



