## General information
- Paper: Altered human oligodendrocyte heterogeneity in multiple sclerosis
- Link: https://www.nature.com/articles/s41586-019-0903-2
- Raw counts downloaded on 2024-01-19
```
wget https://datasets.cellxgene.cziscience.com/57cd8c07-9434-4b4a-923f-8ae275e20b8e.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("57cd8c07-9434-4b4a-923f-8ae275e20b8e.h5ad")
```
- Explore `adata.obs`: 'author_cell_types', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'development_stage'
    - disease: 'normal', 'multiple sclerosis' -> NEED TO subset for normal only
    - organism: Homo sapiens
    - Assay (10x 3' v2)
    - Sex: female and male
    - Tissues: ['brain white matter', 'brain']
    - development_stage: '60-year-old human stage', '35-year-old human stage', '37-year-old human stage', '82-year-old human stage', '64-year-old human stage', '57-year-old human stage', '49-year-old human stage', '44-year-old human stage'
    - How many cell types? 
        - `author_cell_types` (23 cell types): 'COPs', 'Neuron1', 'Neuron2', 'Neuron3', 'Oligo2', 'Astrocytes', 'Oligo6', 'Oligo1', 'Oligo4', 'Oligo3', 'Pericytes', 'Neuron4', 'Neuron5', 'Oligo5', 'ImOlGs', 'Endothelial_cells1', 'Endothelial_cells2', 'OPCs', 'Macrophages', 'Immune_cells', 'Microglia_Macrophages', 'Vasc_smooth_muscle', 'Astrocytes2'
        - `cell_type` (11 cell types): 'oligodendrocyte precursor cell', 'neuron', 'oligodendrocyte', 'astrocyte', 'pericyte', 'endothelial cell', 'progenitor cell', 'macrophage', 'leukocyte', 'microglial cell', 'vascular associated smooth muscle cell'
    - CONCLUSION: keep only `normal`
- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered    feature_name feature_reference feature_biotype feature_length
gene_ids
ENSG00000279457                False          WASH9P    NCBITaxon:9606            gene           1397
ENSG00000228463                False      AP006222.1    NCBITaxon:9606            gene           8224
ENSG00000230021                False  RP11-206L10.17    NCBITaxon:9606            gene           5495
ENSG00000223764                False       LINC02593    NCBITaxon:9606            gene           4152
ENSG00000188976                False           NOC2L    NCBITaxon:9606            gene           5540

adata.var.shape
(20441, 5)
```
- There are 20,441 genes

- Count data is store under the X layer:

```
# the X layer shows raw count
adata.X.A[1:15, 1:15]
array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 2., 1.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 1.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 2., 1.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 1., 0., 0., 2., 0., 0., 1., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]],
      dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                       genes Sample Lesion  ... self_reported_ethnicity        development_stage observation_joinid
Detected                                                    ...
CO28__10X_17_grch38:AAACCTGAGAATCTCCx   2172   CO28   Ctrl  ...                 unknown  60-year-old human stage         ;Ohyox?8=j
CO28__10X_17_grch38:AAATGCCTCACAACGTx   4315   CO28   Ctrl  ...                 unknown  60-year-old human stage         L{bl41M`~m
CO28__10X_17_grch38:AAGGAGCCACGGCGTTx   1998   CO28   Ctrl  ...                 unknown  60-year-old human stage         EPwst(suWV
CO28__10X_17_grch38:ACATACGGTCCAACTAx   1940   CO28   Ctrl  ...                 unknown  60-year-old human stage         fQpg#tSLAm
CO28__10X_17_grch38:ACATCAGAGCCACGCTx   2408   CO28   Ctrl  ...                 unknown  60-year-old human stage         t%h_Q;;Srl

[5 rows x 29 columns]

adata.obs.shape
(17799, 29)
```
- There are 17,799 cells
    - 6,591 are normal

