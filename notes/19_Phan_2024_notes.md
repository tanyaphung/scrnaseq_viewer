## General information
- Paper: Single nuclei transcriptomics in human and non-human primate striatum in opioid use disorder
- Link: https://www.nature.com/articles/s41467-024-45165-7#Abs1
- Raw counts downloaded on 2024-05-30
```
wget https://datasets.cellxgene.cziscience.com/cc8fa27f-1c5a-4ba8-bcf9-17b443fc6aca.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("cc8fa27f-1c5a-4ba8-bcf9-17b443fc6aca.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'author_cell_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 
    - assay: '10x 3' v3'
    - disease: 'normal', 'opiate dependence'
    - organism: 'Homo sapiens'
    - sex: 'female', 'male'
    - tissue: 'caudate nucleus', 'putamen'
    - self_reported_ethnicity: 'European', 'African American'
    - development_stage: '23-year-old human stage', '40-year-old human stage', '41-year-old human stage', '54-year-old human stage', '55-year-old human stage', '35-year-old human stage', '39-year-old human stage', '42-year-old human stage', '37-year-old human stage', '60-year-old human stage'
    - How many cell types?
        - `cell_type` (10 cell type): 'microglial cell', 'oligodendrocyte', 'indirect pathway medium spiny neuron', 'direct pathway medium spiny neuron', 'astrocyte', 'oligodendrocyte precursor cell', 'medium spiny neuron', 'endothelial cell', 'inhibitory interneuron', 'mural cell'
        - `author_cell_type` (12 cell types): "Microglia", "Oligos", "D2-Matrix", "D1-Matrix", "Astrocytes", "Oligos_Pre", "D1-Striosome", "D2-Striosome", "D1/D2-Hybrid", "Endothelial", "Interneuron", "Mural"

- In `adata.var`, gene symbol is under column `feature_name`
```
>>> adata.var.columns
Index(['type', 'feature_is_filtered', 'feature_name', 'feature_reference',
       'feature_biotype', 'feature_length'],
      dtype='object')
>>> adata.var["feature_name"].head()
ENSG00000243485          MIR1302-2HG
ENSG00000186092                OR4F5
ENSG00000238009    ENSG00000238009.6
ENSG00000239906    ENSG00000239906.1
ENSG00000241860    ENSG00000241860.7
Name: feature_name, dtype: category
Categories (31474, object): ['A1BG', 'A1BG-AS1', 'A1CF', 'A2M', ..., 'ZYG11B', 'ZYX', 'ZZEF1', 'ZZZ3']

adata.var.shape
(31474, 6)
```
- There are 31474 genes

- Count data is store under the X layer:

```
adata.X.A[1:10, 1:10]
array([[0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 1.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 2.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                    nCount_RNA  ...  observation_joinid
AAACCCAAGTGCAGGT_1      6472.0  ...          UCXglzX0TE
AAACCCACACTAGGCC_1      4266.0  ...          PuFwsE}7n!
AAACCCACACTGCGTG_1      4813.0  ...          6~v|B2=?#!
AAACCCACAGCGAACA_1     10929.0  ...          -jaTuS0e~z
AAACCCACAGCGGTTC_1      5689.0  ...          rCIes>YJRQ

[5 rows x 58 columns]

adata.obs.shape
(98848, 58)
```
- There are 98848 cells

**IMPORTANT NOTES**
- keep normal only
- separate to caudate nucleus & putamen