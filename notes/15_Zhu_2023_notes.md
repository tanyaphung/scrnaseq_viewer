## General information
- Raw counts downloaded on 2023-12-24 
```
wget https://datasets.cellxgene.cziscience.com/2f17c183-388a-4c08-9adb-a146833e57ab.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("2f17c183-388a-4c08-9adb-a146833e57ab.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'author_cell_type', 'age_group', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity' (unknown), 'development_stage'
    - For developmental stage, there are 6 unique values (column `age_group`): 'adolescence' (group5), 'adulthood' (group6), 'childhood' (group4), 'early fetal' (group1), 'infancy' (group3), 'late fetal' (group2)
        - I will use this to divide the data
        - However, in the column 'development_stage' there are 10 stages and these are a bit more specific: 23rd week post-fertilization human stage, 24th week post-fertilization human stage, 19th week post-fertilization human stage, 18th week post-fertilization human stage, under-1-year-old human stage, 14-year-old human stage, 6-year-old human stage, 4-year-old human stage, 39-year-old human stage, 20-year-old human stage
    - How many cell types? There are 2 colnames for cell type: `cell_type` (13 cell types) and `author_cell_type` (15 cell types)
        - 15 cell types: EN-fetal-late, IN-fetal, VSMC, Endothelial, RG, OPC, IN-CGE, Microglia, IN-MGE, Pericytes, EN-fetal-early, Astrocytes, EN, IPC, Oligodendrocytes
    - For assay, 1 assay: '10x multiome'
    - Only normal
    - For tissue, there are 2 areas: 'cortical plate' and 'dorsolateral prefrontal cortex' but for now I am not separating them

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered   feature_name feature_reference
gene_ids
ENSG00000238009                False   RP11-34P13.7    NCBITaxon:9606  \
ENSG00000241860                False  RP11-34P13.13    NCBITaxon:9606
ENSG00000229905                False  RP11-206L10.4    NCBITaxon:9606
ENSG00000237491                False      LINC01409    NCBITaxon:9606
ENSG00000177757                False         FAM87B    NCBITaxon:9606

                feature_biotype feature_length
gene_ids
ENSG00000238009            gene           3726
ENSG00000241860            gene           7559
ENSG00000229905            gene            456
ENSG00000237491            gene           8413
ENSG00000177757            gene           1947

adata.var.shape
(30113, 5)
```
- There are 30113 genes

- Count data is store under the raw layer:

```
adata.raw.X.A[1:20,1:20]
array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        1., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        1., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        2., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                     author_cell_type   age_group donor_id  nCount_RNA
index
4_AAACAGCCAACACTTG-1    EN-fetal-late  late fetal   LaFet1        3483  \
4_AAACAGCCACCAAAGG-1    EN-fetal-late  late fetal   LaFet1        4863
4_AAACAGCCATAAGTTC-1    EN-fetal-late  late fetal   LaFet1       11069
4_AAACATGCATAGTCAT-1    EN-fetal-late  late fetal   LaFet1        7990
4_AAACATGCATTGTCAG-1    EN-fetal-late  late fetal   LaFet1        6873

adata.obs.shape
(45549, 31)
```
- There are 45549 cells


