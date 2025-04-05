## General information
- Paper: Brain matters: unveiling the distinct contributions of region, age, and sex to glia diversity and CNS function
- Link: https://actaneurocomms.biomedcentral.com/articles/10.1186/s40478-023-01568-z
- Raw counts downloaded on 2024-01-16

```
wget https://datasets.cellxgene.cziscience.com/80a4eb91-dab2-40dc-9856-cb6de4eb954e.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("80a4eb91-dab2-40dc-9856-cb6de4eb954e.h5ad")
```
- Explore `adata.obs`: 'author_cell_type', 'broad_cell_type', 'AgeGroup', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'development_stage'
    - GRCh38; GENCODE 30
    - 2 AgeGroup: 'Old', 'Young'
    - Normal (no disease)
    - Assay (10x 3' v3)
    - Sex: female and male
    - Tissues: ['Brodmann (1909) area 4', 'white matter of cerebellum', 'cervical spinal cord white matter']
    - How many cell types? There are 3 colnames for cell type: 
        - `broad_cell_type` (10 cell types): 'Oligo', 'Neuron_In', 'OPC', 'Neuron_Ex', 'Endothelial-Pericyte', 'Microglia_Macrophages', 'RELN+ neurons', 'Astrocyte', 'Unidentified', 'Neuron'
        - `cell_type` (15 cell types): 'oligodendrocyte', 'GABAergic neuron', 'oligodendrocyte precursor cell', 'glutamatergic neuron', 'endothelial cell of artery', 'capillary endothelial cell', 'mural cell', 'microglial cell', 'cerebellar granule cell', 'vascular associated smooth muscle cell', 'astrocyte', 'leukocyte', 'differentiation-committed oligodendrocyte precursor', 'central nervous system macrophage', 'neuron'
        - `author_cell_type` (60 cell type)
    - CONCLUSION: separate out by `tissue` and `AgeGroup`
- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered  vst.mean  vst.variance  vst.variance.expected  ...   feature_name feature_reference feature_biotype feature_length
ENSG00000237613                 True       NaN           NaN                    NaN  ...        FAM138A    NCBITaxon:9606            gene           1219
ENSG00000238009                False  0.004792      0.005593                0.00659  ...   RP11-34P13.7    NCBITaxon:9606            gene           3726
ENSG00000239945                 True       NaN           NaN                    NaN  ...   RP11-34P13.8    NCBITaxon:9606            gene           1319
ENSG00000239906                 True       NaN           NaN                    NaN  ...  RP11-34P13.14    NCBITaxon:9606            gene            323
ENSG00000284733                 True       NaN           NaN                    NaN  ...         OR4F29    NCBITaxon:9606            gene            939

adata.var.shape
(33784, 9)
```
- There are 33,784 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata.X.A[1:25, 1:25]
array([[0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.5901919 ,
        0.        , 0.        , 0.        , 1.0078369 ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 2.6348069 ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 2.0779276 ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 1.5592016 ],
# the raw count is in the raw.X layer
adata.raw.X.A[1:20, 1:20]
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
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 2., 0., 0., 0., 0.,
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
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                                 mapped_reference_assembly mapped_reference_annotation  ...        development_stage observation_joinid
BA4_possorted_genome_bam_A624Z:AAAGTGAGTCCACGCAx                    GRCh38                  GENCODE 30  ...  63-year-old human stage         ARTth)}Z=m
BA4_possorted_genome_bam_A624Z:ACATTTCGTTAAGTCCx                    GRCh38                  GENCODE 30  ...  63-year-old human stage         TfWG6eLEV+
BA4_possorted_genome_bam_A624Z:AACCAACCATCACGGCx                    GRCh38                  GENCODE 30  ...  63-year-old human stage         0m+pbH~1@j
BA4_possorted_genome_bam_A624Z:ACAGCCGAGTGCACTTx                    GRCh38                  GENCODE 30  ...  63-year-old human stage         mPafjB5z$5
BA4_possorted_genome_bam_A624Z:AACCCAAGTGTTTGCAx                    GRCh38                  GENCODE 30  ...  63-year-old human stage         !Lw?i1p6jS

adata.obs.shape
(45528, 49)
```
- There are 45,528 cells

