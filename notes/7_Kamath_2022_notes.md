## General information
- Paper: Single-cell genomic profiling of human dopamine neurons identifies a population that selectively degenerates in Parkinson's disease
- Link: https://pubmed.ncbi.nlm.nih.gov/35513515/
- Raw counts downloaded on 2024-05-30

```
wget https://datasets.cellxgene.cziscience.com/ba41e78a-6b71-40d6-9b40-50a61800cb5c.h5ad
mv ba41e78a-6b71-40d6-9b40-50a61800cb5c.h5ad oligodendrocytes.h5ad

wget https://datasets.cellxgene.cziscience.com/51affad7-08ed-43b3-93ad-56640f4b8910.h5ad
mv 51affad7-08ed-43b3-93ad-56640f4b8910.h5ad nurr_negative.h5ad

wget https://datasets.cellxgene.cziscience.com/a2dc03fb-8820-4d3d-a7ac-ab40686bab3c.h5ad
mv a2dc03fb-8820-4d3d-a7ac-ab40686bab3c.h5ad non_DA_neurons.h5ad

wget https://datasets.cellxgene.cziscience.com/8cfb7c71-ba45-4b1e-8df9-c84966f069ee.h5ad
mv 8cfb7c71-ba45-4b1e-8df9-c84966f069ee.h5ad nurr_positive.h5ad

wget https://datasets.cellxgene.cziscience.com/a0559ab2-c0d0-41e0-af27-a3ff87d3ae28.h5ad
mv a0559ab2-c0d0-41e0-af27-a3ff87d3ae28.h5ad astrocytes.h5ad

wget https://datasets.cellxgene.cziscience.com/c480e527-5725-4699-bd8a-e09535b23ba8.h5ad
mv c480e527-5725-4699-bd8a-e09535b23ba8.h5ad microglia.h5ad

wget https://datasets.cellxgene.cziscience.com/e6b3b510-a431-4626-a9ac-76e7b806aaa2.h5ad
mv e6b3b510-a431-4626-a9ac-76e7b806aaa2.h5ad DA_neurons.h5ad

wget https://datasets.cellxgene.cziscience.com/0566a345-7bde-4107-9976-3f7976c0c5ca.h5ad
mv 0566a345-7bde-4107-9976-3f7976c0c5ca.h5ad endothelial.h5ad

wget https://datasets.cellxgene.cziscience.com/675b276f-6b88-4c65-aa37-7d1df922b705.h5ad #opc
mv 675b276f-6b88-4c65-aa37-7d1df922b705.h5ad opc.h5ad
```

## Data exploration
### Example for OPC
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("opc.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'author_cell_type', 'cell_type', 'assay', 'disease',
       'organism', 'sex', 'tissue', 'self_reported_ethnicity',
       'development_stage'
    - assay: '10x 3' v3'
    - disease: 'normal', 'Lewy body dementia', 'Parkinson disease'
    - organism: 'Homo sapiens'
    - sex: 'male', 'female'
    - tissue: 'substantia nigra pars compacta'
    - self_reported_ethnicity: 'unknown'
    - development_stage: '79-year-old human stage', '90-year-old human stage', '82-year-old human stage', '91-year-old human stage', '92-year-old human stage', ..., '74-year-old human stage', '75-year-old human stage', '76-year-old human stage', '83-year-old human stage', '78-year-old human stage'
    - How many cell types?
        - `cell_type` (1 cell type): 'oligodendrocyte precursor cell'
        - `author_cell_type` (5 cell types): 'OPC_CACNG4', 'OPC_HOXD3', 'OPC_ADM', 'OPC_KIAA0040', 'OPC_MDFI'

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered       feature_name  ... feature_biotype feature_length
ENSG00000237613                False            FAM138A  ...            gene           1219
ENSG00000186092                False              OR4F5  ...            gene           2618
ENSG00000238009                False  ENSG00000238009.6  ...            gene           3726
ENSG00000239945                False  ENSG00000239945.1  ...            gene           1319
ENSG00000239906                False  ENSG00000239906.1  ...            gene            323

[5 rows x 5 columns]

adata.var.shape
(34010, 5)
```
- There are 34010 genes

- Count data is store under the X layer:

```
adata.X.A[1:10, 1:10]
array([[0., 1., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 1.],
       [0., 0., 0., 0., 0., 0., 0., 0., 3.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                               author_cell_type  ... observation_joinid
index                                                            ...
pPDsHSrSNxi3298d200429DAPIA_AACCAACTCGTTCCTG-1       OPC_CACNG4  ...         BZ{mr<TG5y
pPDsHSrSNxi3298d200429DAPIA_TCCTTTCCAGCCGTCA-1       OPC_CACNG4  ...         cm}xmhB~8L
pPDsHSrSNxi3298d200429DAPIA_ACTATCTTCAGACCGC-1       OPC_CACNG4  ...         ;SFZzSbjwL
pPDsHSrSNxi3298d200429DAPIA_TCAGTCCTCAAACCCA-1       OPC_CACNG4  ...         _8<j>o>fkF
pPDsHSrSNxi3298d200429DAPIA_CACAACACAAGACCTT-1       OPC_CACNG4  ...         +nR5^kMXP`

[5 rows x 25 columns]

adata.obs.shape
(13691, 25)
```
- There are 13691 cells

**IMPORTANT NOTES**
- keep normal only
- For this dataset, each h5ad file that I downloaded is for one cell type. Therefore, I will subset each cell type file first for normal and then concatenate