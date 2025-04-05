## General information
- Paper: Spatial and cell type transcriptional landscape of human cerebellar development
- Link: https://www.nature.com/articles/s41593-021-00872-y#Abs1
- Raw counts downloaded on 2024-01-07
```
wget https://datasets.cellxgene.cziscience.com/14c2c322-e424-47b3-a7af-dd19df10b644.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("14c2c322-e424-47b3-a7af-dd19df10b644.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'figure_clusters', 'experiment', 'fig_cell_type',
       'donor_id', 'cell_type', 'assay', 'disease', 'sex',
       'tissue', 'development_stage'
    - Normal (no disease)
    - Assay (SPLiT-seq only)
    - Sex: female and male
    - 13 donors
    - Development_stage: 9 total (however, from Fig1 of the paper it indicates 10): 5003+2329+20364+7747+11195+15549+1601+5177+143
        - 9th week post-fertilization human stage (5003 cells)
        - 10th week post-fertilization human stage (2329 cells)
        - 11th week post-fertilization human stage (20364 cells)
        - 12th week post-fertilization human stage (7747 cells)
        - 14th week post-fertilization human stage (11195 cells)
        - 17th week post-fertilization human stage (15549 cells)
        - 18th week post-fertilization human stage (1601 cells)
        - 20th week post-fertilization human stage (5177 cells)
        - 21st week post-fertilization human stage (143 cells)
        - (missing 16 week)
    - How many cell types? There is 2 colname for cell type: 
        - `cell_type` (18 cell types): glial cell, cerebellar granule cell precursor, Purkinje cell, meningeal macrophage, CNS interneuron, Bergmann glial cell, astrocyte, oligodendrocyte precursor cell, ependymal cell, granule cell, interneuron, microglial cell, inhibitory interneuron, choroid plexus epithelial cell, GABAergic interneuron, brainstem motor neuron, endothelial cell, pericyte
        - `fig_cell_type` (21 cell types): H-Glia, H-RL, H-GCP, H-PC, H-Meninges, H-eCN/UBC, H-BG, H-Ast, H-OPC, H-Ast/Ependymal, H-GN, H-PIP, H-Microglia, H-iCN, H-Choroid, H-MLI, H-Brainstem, H-BS Choroid/Ependymal, H-Endothelial, H-Pericytes, H-Committed OPC
    - 3 experiments.  

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                     name          ensg_id  feature_is_filtered
ensg_id
ENSG00000108984    MAP2K6  ENSG00000108984                False  \
ENSG00000159399       HK2  ENSG00000159399                False
ENSG00000033867    SLC4A7  ENSG00000033867                False
ENSG00000254093     PINX1  ENSG00000254093                False
ENSG00000169905  TOR1AIP2  ENSG00000169905                False

                          feature_name feature_reference feature_biotype
ensg_id
ENSG00000108984                 MAP2K6    NCBITaxon:9606            gene  \
ENSG00000159399                    HK2    NCBITaxon:9606            gene
ENSG00000033867                 SLC4A7    NCBITaxon:9606            gene
ENSG00000254093  PINX1_ENSG00000254093    NCBITaxon:9606            gene
ENSG00000169905               TOR1AIP2    NCBITaxon:9606            gene

                feature_length
ensg_id
ENSG00000108984          16660
ENSG00000159399           5918
ENSG00000033867           9616
ENSG00000254093           2804
ENSG00000169905          18876

adata.var.shape
(26361, 7)
```
- There are 26,361 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata.X.A[1:10, 1:10]
array([[0.       , 0.       , 0.       , 0.       , 0.       , 2.5289779,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 2.0515127, 0.       , 3.0605357,
        0.       , 0.       , 0.       ],
       [0.       , 2.429391 , 0.       , 0.       , 0.       , 4.6493754,
        0.       , 2.429391 , 0.       ],
       [0.       , 2.185927 , 0.       , 0.       , 0.       , 2.8212473,
        0.       , 0.       , 0.       ],
       [0.       , 3.1819913, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 3.0579545,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 4.2292476,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 1.9571365,
        0.       , 0.       , 0.       ]], dtype=float32)
# the raw count is in the raw.X layer
adata.raw.X.A[1:10, 1:10]
array([[ 0.,  0.,  0.,  0.,  0.,  2.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  1.,  0.,  3.,  0.,  0.,  0.],
       [ 0.,  1.,  0.,  0.,  0., 10.,  0.,  1.,  0.],
       [ 0.,  1.,  0.,  0.,  0.,  2.,  0.,  0.,  0.],
       [ 0.,  2.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  2.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0., 14.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  1.,  0.,  0.,  0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
       orig.ident  nCount_RNA  nFeature_RNA  percent.mito   S.Score
V3_1          01k       617.0           482      0.012966 -0.012428  \
V224_1        01k      1733.0          1147      0.010964 -0.034455
V364_1        01k      1475.0           974      0.027797 -0.026429
V372_1        01k       966.0           718      0.010352  0.060085
V434_1        01k      1266.0           977      0.030806 -0.029667

adata.obs.shape
(69174, 35)
```
- There are 69,174 cells
- CONCLUSION: will separate out into 9 developmental stage and will not separate out by experiment

