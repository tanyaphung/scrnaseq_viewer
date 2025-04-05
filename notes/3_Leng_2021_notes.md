## General information
- Paper: Molecular characterization of selectively vulnerable neurons in Alzheimer's Disease
- Link: https://www.nature.com/articles/s41593-020-00764-7
- Raw counts downloaded on 2024-05-30

```
wget https://datasets.cellxgene.cziscience.com/a6a7e2af-d0a2-4091-beb5-f65b2b9d59a1.h5ad #superior frontal gyrus
mv a6a7e2af-d0a2-4091-beb5-f65b2b9d59a1.h5ad sfg.h5ad

wget https://datasets.cellxgene.cziscience.com/44703324-8d6c-44c2-958d-8d16f5eb4054.h5ad #entorhinal cortex
mv 44703324-8d6c-44c2-958d-8d16f5eb4054.h5ad ec.h5ad
```

## Data exploration
### superior frontal gyrus
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("sfg.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'clusterAssignment',
       'clusterCellType', 'cell_type', 'assay', 'disease', 'organism', 'sex',
       'tissue', 'self_reported_ethnicity', 'development_stage'
    - clusterAssignment: SFG:Oligo.2, SFG:Micro, SFG:Oligo.1, SFG:OPC, SFG:Exc.1, SFG:Astro.2, SFG:Exc.2, SFG:Inh.4, SFG:Exc.3, SFG:Astro.1, SFG:Exc.4, SFG:Exc.6, SFG:Inh.3, SFG:Exc.5, SFG:Inh.2, SFG:Inh.1, SFG:Endo, SFG:Exc.7
    - clusterCellType: 'Oligo', 'Micro', 'OPC', 'Exc', 'Astro', 'Inh', 'Endo'
    - assay: '10x 3' v2'
    - disease: 'normal', 'Alzheimer disease'
    - organism: 'Homo sapiens'
    - sex: 'male'
    - tissue: 'superior frontal gyrus'
    - self_reported_ethnicity: unknown
    - development_stage: '60-year-old human stage', '50-year-old human stage', '71-year-old human stage', '72-year-old human stage', '87-year-old human stage', '80 year-old and over human stage', '77-year-old human stage', '82-year-old human stage'
    - How many cell types?
        - `cell_type` (7 cell types): 'oligodendrocyte', 'mature microglial cell', 'oligodendrocyte precursor cell', 'glutamatergic neuron', 'mature astrocyte', 'GABAergic neuron', 'endothelial cell'

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered       feature_name  ... feature_biotype feature_length
ENSG00000100568                False              VTI1B  ...            gene           5768
ENSG00000101440                False               ASIP  ...            gene            845
ENSG00000249847                False  ENSG00000249847.1  ...            gene            646
ENSG00000136630                False                HLX  ...            gene           5629
ENSG00000231731                False  ENSG00000231731.7  ...            gene           5510

[5 rows x 5 columns]

adata.var.shape
(32743, 5)
```
- There are 32,743 genes

- Count data is store under the X layer:

```
adata.X.A[1:10, 1:10]
array([[0., 0., 0., 0., 0., 0., 1., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 1., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 1., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                      SampleID donor_id  ...        development_stage observation_joinid
SFG2_AAACCTGAGATGGCGT     SFG2        2  ...  60-year-old human stage         Rw$9Oe)Es<
SFG2_AAACCTGAGCGATCCC     SFG2        2  ...  60-year-old human stage         Agc~%jD4^-
SFG2_AAACCTGAGGAATCGC     SFG2        2  ...  60-year-old human stage         ILz&V-ZSaz
SFG2_AAACCTGAGGATGCGT     SFG2        2  ...  60-year-old human stage         5wsZB%4>$R
SFG2_AAACCTGAGGCACATG     SFG2        2  ...  60-year-old human stage         -j4TM>0zT?

[5 rows x 30 columns]

adata.obs.shape
(63608, 30)
```
- There are 63608 cells

**IMPORTANT NOTES**
- keep normal only

### entorhinal cortex
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("ec.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'clusterAssignment',
       'clusterCellType', 'cell_type', 'assay', 'disease', 'organism', 'sex',
       'tissue', 'self_reported_ethnicity', 'development_stage'
    - clusterAssignment: EC:Exc.5, EC:Exc.3, EC:Exc.1, EC:OPC, EC:Inh.1, EC:Exc.2, EC:Micro, EC:Astro, EC:Inh.2, EC:Oligo, EC:Exc.4, EC:Inh.3, EC:Endo
    - clusterCellType: 'Exc', 'OPC', 'Inh', 'Micro', 'Astro', 'Oligo', 'Endo'
    - assay: '10x 3' v2'
    - disease: 'normal', 'Alzheimer disease'
    - organism: 'Homo sapiens'
    - sex: 'male'
    - tissue: 'entorhinal cortex'
    - self_reported_ethnicity: unknown
    - development_stage: '60-year-old human stage', '50-year-old human stage', '71-year-old human stage', '72-year-old human stage', '87-year-old human stage', '80 year-old and over human stage', '77-year-old human stage', '82-year-old human stage'
    - How many cell types?
        - `cell_type` (7 cell types): 'glutamatergic neuron', 'oligodendrocyte precursor cell', 'GABAergic neuron', 'mature microglial cell', 'mature astrocyte', 'oligodendrocyte', 'endothelial cell'

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered       feature_name  ... feature_biotype feature_length
ENSG00000100568                False              VTI1B  ...            gene           5768
ENSG00000101440                False               ASIP  ...            gene            845
ENSG00000249847                False  ENSG00000249847.1  ...            gene            646
ENSG00000136630                False                HLX  ...            gene           5629
ENSG00000231731                False  ENSG00000231731.7  ...            gene           5510

[5 rows x 5 columns]

adata.var.shape
(32743, 5)
```
- There are 32,743 genes

- Count data is store under the X layer:

```
adata.X.A[1:10, 1:10]
adata.X.A[1:10, 1:10]
array([[0., 0., 0., 0., 0., 0., 1., 0., 0.],
       [0., 0., 0., 0., 1., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [1., 0., 0., 0., 0., 0., 2., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                     SampleID donor_id  ...        development_stage observation_joinid
EC2_AAACCTGAGGATGCGT      EC2        2  ...  60-year-old human stage         SVJIToe1C=
EC2_AAACCTGAGTCAATAG      EC2        2  ...  60-year-old human stage         VWzMofWup%
EC2_AAACCTGCATACTCTT      EC2        2  ...  60-year-old human stage         ZKMx*@_bx0
EC2_AAACCTGGTCCAACTA      EC2        2  ...  60-year-old human stage         -GV819~6Fr
EC2_AAACCTGGTTATCACG      EC2        2  ...  60-year-old human stage         oAFf|pFl@q

[5 rows x 30 columns]

adata.obs.shape
(42528, 30)
```
- There are 42528 cells

**IMPORTANT NOTES**
- keep normal only