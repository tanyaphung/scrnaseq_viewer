# scrnaseq_viewer
This repository provides an overview of single cell RNAseq data that I have processed

# Tissue: Brain
## Developmental stage: adult
1. Jorstad et al. 2023
- Title: Comparative transcriptomics reveals human-specific cortical features
- Link: https://www.science.org/doi/10.1126/science.ade9516
- Specific region: middle temporal gyrus (MTG)
- Information from the observation dataframe:
  - Cell type annotations are stored in 3 variables: Cluster (151 annotations), Subclass (24 annotations), and cell_type (18)
  - Developmental stage: all adults (human adult stage, 29-year-old human stage, 42-year-old human stage, 43-year-old human stage, 50-year-old human stage, and 60-year-old human stage)
  - Ethnicity is unknown
  - Total of 7 donors
  - 2 assays: Smart-seq v4 and 10x 3' v3
  - Normal (not from disease)
  - Both male (2) and female (5)
- Cell type labels are similar to Allen Brain Atlas
  - Inhibitory neurons
    - 5 CGE-derived subclasses (LAMP5 LHX6 LAMP5, VIP, PAX6, and SNCG) expressing the marker ADARB2
    - 4 MGE-derived subclasses (Chandelier, PVALB, SST, and SST CHODL) expressing LHX6
  - Excitatory neurons
    - 5 IT-projecting subclasses (L2/3 IT, L4 IT, L5 IT, L6 IT, and L5/6 IT CAR3)
    - 4 deep-layer non-IT-projecting subclasses (L5 ET, L5/6 NP, L6b, and L6 CT)
  - Non-neuronal cells: astrocytes, OPCs, oligodendrocytes, microglia/PVMs, endothelia cells, and VLMCs
