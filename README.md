![image](https://github.com/tanyaphung/scrnaseq_viewer/assets/10180091/a66fbb1f-2eee-439e-a1c0-e021e855179b)# scrnaseq_viewer
This repository provides an overview of single cell RNAseq data that I have processed

# Tissue: Brain
## Developmental stage: adult
1. Jorstad et al. 2023
- Title: Comparative transcriptomics reveals human-specific cortical features
- Link: https://www.science.org/doi/10.1126/science.ade9516
- Short summary: Jorstad et al. 2023 sequenced the MTG region of the adult brain which resulted in about 150,000 cells total
- Specific regions: middle temporal gyrus (MTG)
- Total number of cells: 156285
- Information from the observation dataframe:
  - Cell type annotations are stored in 3 variables: Cluster (151 annotations), Subclass (24 annotations), and cell_type (18)
  - Developmental stage: all adults (human adult stage, 29-year-old human stage, 42-year-old human stage, 43-year-old human stage, 50-year-old human stage, and 60-year-old human stage)
  - Ethnicity is unknown
  - Total of 7 donors
  - 2 assays: Smart-seq v4 (14503 cells) and 10x 3' v3 (141782 cells)
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

## Developmental stage: multiple stages (prenatal, postnatal, and adult)
1. Velmeshev et al. 2023
- Title: Single-cell analysis of prenatal and postnatal human cortical development
- Link: https://www.science.org/doi/10.1126/science.adf0834
- Short summary:
- Specific regions: ganglionic eminences (major source of cortical interneurons) and from the cortex (prefrontal, cingulate, temporal, insular, and motor cortices)
- Total number of cells: in this study, there were a total of 413,682 nuclei. After excluding some specific cluster after QC, there were 358663 nuclei left at the end. Then, the authors integrate with existing datasets which results in 709,372 nuclei.
- Information from the observation dataframe:
  - Data from this study (Velmeshev) was combined with 3 other studies: Ramos, Herring, and Trevino
  - Assay: all 10x but different flavors: 10x 3' v2, 10x 3' v3, 10x multiome
  - Normal (not from disease)
  - Both male and female
  - Developmental stage: there are 40 unique values but from Figure 1 of the paper, they could be grouped in 8 groups: 2nd trimester, 3rd trimester, 0-1 years, 1-2 years, 2-4 years, 4-10 years, 10-20 years, and Adult. 
- Cell type labels:
  - 6: native cell, astrocyte, oligodendrocyte, microglial cell, neural cell, oligodendrocyte precursor cell
