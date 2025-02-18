# scrnaseq_viewer
This repository provides an overview of single cell RNAseq data that I have processed as part of Phung et al. (in preparation).

# Setting up conda environment
- yaml file is provided in `resources/scrna_env.yml`
```
conda env create --name scrna --file=resources/scrna_env.yml
```

##### Table of Contents  
[Tissue: Brain](#tissue-brain)  

  [Developmental stage: adult](#developmental-stage-adult)
  
  [Developmental stage: prenatal and early postnatal](#developmental-stage-prenatal-and-early-postnatal)
  
  [Developmental stage: child](#developmental-stage-child)
  
  [Developmental stage: multiple stages (prenatal, postnatal, and adult)](#developmental-stage-multiple-stages-prenatal-postnatal-and-adult)


# Tissue: Brain
## Developmental stage: adult
(Note: the order of the datasets described here follows the ordering from the manuscript in preparation - more details to follow)
#### 13. Jorstad et al. 2023 MTG
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

#### 16. Jorstad et al. 2023 Neocortex
- Title: Transcriptomic cytoarchitecture reveals principles of human neocortex organization
- Link: [https://www.science.org/doi/10.1126/science.ade9516](https://pubmed.ncbi.nlm.nih.gov/37824655/)
- Short summary: To better understand cortical areal specilization, Jorstad et al. 2023 sequenced the neocortex of the adult brain.
- 8 areas spanning cortical structural variation: primary motor cortex (M1), primary somatosensory cortex (S1), primary auditory cortex (A1), primary visual cortex (V1), dorsolateral prefrontal cortex (DFC), anterior cingulate cortex (ACC), middle temporal gyrus (MTG), and angular gyrus (AnG)
- Total number of cells: ~1.15 million from 10x and ~49k from Smart-seq
 
#### 6. Bakken et al. 2021
**TL;DR:** primary motor cortex, adult
- Title: Comparative cellular analysis of motor cortex in human, marmoset and mouse
- Link: https://www.nature.com/articles/s41586-021-03465-8#Abs1
- Short summary: Bakken et al. 2021 sequenced the primary motor cortex region of the adult brain from different species (humans, marmoset monkey, and mice)
- Specific region: primary motor cortex M1
- Total number of cells: ~20k for humans only
- Assay: 10x 3' v3

#### 14. Siletti et al. 2023
- Title: Transcriptomic diversity of cell types across the adult human brain
- Link: https://pubmed.ncbi.nlm.nih.gov/37824663/
- Short summary: 
  - Isolated postmortem tissue from 3 donors
  - Enriched for neurons from ~100 locations across the forebrain (cerebral cortex, hippocampus, cerebral nuclei, hypothalamus, and thalamus), midbrain, and hindbrain (pons, medulla, and cerebellum)
  - ~3M cells  31 superclusters, 461 clusters, and 3313 subclusters
- Preprocessing scRNA steps
  - `code/preprocessing/14_Siletti_AdultHumanBrain/qc_scrna_siletti_adulthumanbrain.py`
    - I used a snakemake pipeline to run this python script in parallel: `code/preprocessing/14_Siletti_AdultHumanBrain/snakefile_qc_scrna_siletti.smk`
  - update cell type columns (this is because when this data was originally processed, we did not come to the conclusion in our group that we would rename the first level of cell type annotation to cell_type_level_1)
    - `code/preprocessing/14_Siletti_AdultHumanBrain/update_celltype.py` (run with `code/preprocessing/14_Siletti_AdultHumanBrain/run_update_celltype.sh`)

## Developmental stage: prenatal and early postnatal
#### 1. Smith et al. 2021 
**TL;DR:** neocortex; midgestational and infant
- Title: Early role for a Na+,K+-ATPase (ATP1A3) in brain development
- Link: https://pubmed.ncbi.nlm.nih.gov/34161264/
- Short summary: In this study Smith et al. 2021 is mostly interested in the gene ATP1A3 because variants in this gene could be linked to polymicrogyria, a developmental malformation of the cerebral cortex characterized by abnormal folding and laminar organization. 
- Specific region: neocortex
- Midgestational:
  - Total number of human normal cells: 118,647
  - Assay: Dropseq
  - Only normal
  - Tissues: 'parietal lobe', 'hippocampal formation', 'primary visual cortex', 'medial ganglionic eminence', 'caudal ganglionic eminence', 'orbitofrontal cortex', 'anterior cingulate cortex'
  - Cell type lables: 'glial cell', 'GABAergic neuron', 'glutamatergic neuron', 'neural progenitor cell'
- Infant:
  - Total number of human normal cells: 51,878
  - Assay: Dropseq
  - Only normal
  - Tissues: 'prefrontal cortex', 'temporal lobe', 'parietal lobe', 'primary visual cortex'
  - Cell type lables: 'glial cell', 'GABAergic neuron', 'glutamatergic neuron'
 
#### 4. Aldinger et al. 2021
**TL;DR:** cerebellum; 9-21 wpc
- Title: Spatial and cell type transcriptional landscape of human cerebellar development
- Link: https://pubmed.ncbi.nlm.nih.gov/34140698/
- Short summary: Since the cerebellum is not well represented in previous bulk and single-cell transcriptomic studies of the developing human brain, the developmental origins of the cerebellar neuroanatomy is not well understood. Therefore, the goal of this study is to perform both bulk RNAseq and scRNAseq across 9 developmental stages prenatally ranging from 9 weeks to 21 weeks post fertilization
- Total number of human normal cells: 69174
- Assay: SPLiT-seq
- Tissue: cerebella

#### 5. Bhaduri et al. 2021
**TL;DR:** 6 neocortical areas; 2nd trimester
- Title: An atlas of cortical arealization identifies dynamic molecular signatures
- Link: [https://pubmed.ncbi.nlm.nih.gov/34140698/](https://pubmed.ncbi.nlm.nih.gov/34616070/)
- Short summary: How early patterns result in finer and more discrete spatial differences remains poorly understood. To contribute to better understand this, the authors sequenced 10 major brain structures including 6 neocortical areas during the second trimester
- Total number of human normal cells: 698,820 but only 457,965 cells are available from the download
- Assay: 10x 3' v2
- Tissue: neocortex
 
## Developmental stage: child
#### 21. Torres-Flores
**TL;DR:** cerebellum; child_age
- Short summary: 
- Total number of human normal cells: 
- Assay: 10x 3' v3
- Tissue: cerebellum

## Developmental stage: multiple stages (prenatal, postnatal, and adult)
#### 17. Velmeshev et al. 2023
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

#### 18. Sepp et al. 2023
- Title: Cellular development and evolution of the mammalian cerebellum
- Link: https://www.nature.com/articles/s41586-023-06884-x#Abs1
- Short summary: Sepp et al. 2023 sequenced the cerebellum in mouse, human, and opossum at different developmental stages. There were ~400k cells in total but only ~160k normal human cells
- Specific region: cerebellum
- Total number of human normal cells: 163,283
- Information from the observation dataframe:
  - Assay: there seems to be 2 versions of 10x: 10x 3' v2, 10x 3' v3. For now I won't separate them but keep in mind there might be some differences. See Extended Data Fig. 1.
  - Data included diseases and normal so I only keep the cells coming from normal (healthy) samples
  - Both male and female
  - Developmental stage: from the column author_stage, there are 18 unique values but I am going to group them into 11 groups: CS18, CS19, CS22, 9 wpc, 11 wpc, 17 wpc, 20 wpc, newborn, infant, toddler, and adult
  - Cell type lables:
    - 25: native cell, cerebellar granule cell, glutamatergic neuron, macroglial cell, Purkinje cell neuroblast, GABAergic neuron, cerebellar granule cell precursor, interneuron, unipolar brush cell, microglial cell, progenitor cell, gllioblast, noradrenergic cell, central nervous system macrophage, brain vascular cell, erythroid lineage cell, leukocyte, oligodendrocyte precursor cell, T cell, meningeal macrophage, Bergmann glial cell, immature astrocyte, oligodendrocyte

#### 15. Zhu et al. 2023
- Title: Multi-omic profiling of the developing human cerebral cortex at the single-cell level
- Link: https://www.science.org/doi/10.1126/sciadv.adg3754
- Short summary: Zhu et al. 2023 sequenced the neocortex using both snRNAseq and scATACseq techonologies (simutaneous multi-omics single-cell profiling in the developing brain)
- Specific region: neocortex
- Total number of human normal cells: 45,549
- Information from the observation dataframe:
  - Assay: 10x multiome
  - Only normal
  - Developmental stage (6 stages): 18-19 GW, 23-24 GW, 0 year, 4-6 years, 14 years, 20-39 years
  - Cell type lables:
    - 13 cell types from the column `cell_type`
    - 15 cell types from the column `author_cell_type`: EN-fetal-late, IN-fetal, VSMC, Endothelial, RG, OPC, IN-CGE, Microglia, IN-MGE, Pericytes, EN-fetal-early, Astrocytes, EN, IPC, Oligodendrocytes


