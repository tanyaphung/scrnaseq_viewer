# Documentation
## General information
- Paper: no paper linked to this dataset
- Raw counts downloaded on 2024-01-10 by Tanya Phung
```
pwd; date
/home/tphung/tphung_proj/Preprocessing_scRNA/data/4_TorresFlores_Cerebellum_Human
Wed Jan 10 13:29:56 CET 2024
wget https://datasets.cellxgene.cziscience.com/7dd2aaf8-1714-49bc-a434-671e8465ace0.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("7dd2aaf8-1714-49bc-a434-671e8465ace0.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'customclassif', 'cell_type', 'assay',
       'disease', 'organism', 'sex', 'tissue', 'development_stage'
    - Disease: 'normal', 'pilocytic astrocytoma' -> need to subset for normal only
    - Assay (10x 3' v3)
    - Sex: male and female
    - Organism: Homo sapiens
    - Tissue: cerebellum
    - Development stage: 'child stage'
    - How many cell types? 
        - `cell_type` (9 cell types): 'neuron', 'oligodendrocyte', 'macrophage', 'oligodendrocyte precursor cell', 'endothelial cell', 'Bergmann glial cell', 'astrocyte', 'interneuron', 'T cell'
        - `customclassif` (10 cell types): 'Inhibitory neurons', 'Excitatory neurons', 'Oligodendrocytes', 'Macrophage_CCL4 high', 'Oligodendrocyte precursor cells', 'Endothelial cells', 'Astrocyte (Bergmann glia)', 'Astrocytes', 'Interneurons', 'T cells' -> use this as level 1

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 vst.mean  vst.variance  vst.variance.expected
gene_ids
ENSG00000237491  0.025956      0.033645               0.034847  \
ENSG00000228794  0.014086      0.015684               0.018235
ENSG00000188976  0.020877      0.022743               0.027514
ENSG00000188290  0.021495      0.031584               0.028389
ENSG00000187608  0.037321      0.042888               0.052613

                 vst.variance.standardized  vst.variable  feature_is_filtered
gene_ids
ENSG00000237491                   0.965524         False                False  \
ENSG00000228794                   0.860144         False                False
ENSG00000188976                   0.826598         False                False
ENSG00000188290                   1.112562         False                False
ENSG00000187608                   0.815171         False                False

                feature_name feature_reference feature_biotype feature_length
gene_ids
ENSG00000237491    LINC01409    NCBITaxon:9606            gene           8413
ENSG00000228794    LINC01128    NCBITaxon:9606            gene          15682
ENSG00000188976        NOC2L    NCBITaxon:9606            gene           5540
ENSG00000188290         HES4    NCBITaxon:9606            gene           1118
ENSG00000187608        ISG15    NCBITaxon:9606            gene            867

adata.var.shape
(28801, 10)
```
- There are 28,801 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata.X.A[1:10, 1:10]
array([[0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        2.7123542, 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        2.563105 , 0.       , 0.       ]], dtype=float32)
# the raw count is in the raw.X layer
adata.raw.X.A[1:10, 1:10]
array([[0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 1., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 1., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                       seq_folder    nUMI  nGene          astro  percent.mt
index
1.1_AAACGAACAACGACAG-1       m184  2306.0   1358  Astrocytoma 1    2.601908  \
1.1_AAACGAACATTGAAGA-1       m184   711.0    552  Astrocytoma 1    1.265823
1.1_AAACGCTCACCAGTAT-1       m184   820.0    593  Astrocytoma 1    1.219512
1.1_AAAGGGCGTTACAGCT-1       m184   525.0    432  Astrocytoma 1    2.476190
1.1_AAAGGTACAACTCGAT-1       m184  1021.0    722  Astrocytoma 1    1.077375

adata.obs.shape
(35637, 35)
```
- There are 35,637 cells

- Conclusion: just need to subset for normal samples