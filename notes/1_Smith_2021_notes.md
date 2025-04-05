## General information
- Paper: https://www.pnas.org/doi/full/10.1073/pnas.2023333118#supplementary-materials
- Link: https://cellxgene.cziscience.com/collections/e02201d7-f49f-401f-baf0-1eb1406546c0
- Raw counts downloaded on 2023-12-25
```
# Midgestational human neocortex
wget https://datasets.cellxgene.cziscience.com/f98941c2-9747-49e5-8b8a-38b183658d85.h5ad

# Infant human neocortex
wget https://datasets.cellxgene.cziscience.com/30b1c353-d4e1-4f34-9577-54dcc19ab4ff.h5ad
```

## Data exploration
### Midgestational human neocortex
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("f98941c2-9747-49e5-8b8a-38b183658d85.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'cell_type', 'assay', 'disease', 'sex', 'tissue', 'development_stage'
    - Normal (no disease)
    - Assay (Drop-seq only)
    - Sex (unknown)
    - Development_stage: 22nd week post-fertilization human stage
    - How many cell types? There is 1 colname for cell type: `cell_type` (4 cell types): 'glial cell', 'GABAergic neuron', 'glutamatergic neuron', 'neural progenitor cell'
    - For tissue, there are 7 areas: 'parietal lobe', 'hippocampal formation', 'primary visual cortex', 'medial ganglionic eminence', 'caudal ganglionic eminence', 'orbitofrontal cortex', 'anterior cingulate cortex'

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                            name  feature_is_filtered   feature_name
index
ENSG00000241860  ENSG00000241860                False  RP11-34P13.13  \
ENSG00000228463  ENSG00000228463                False     AP006222.1
ENSG00000228327  ENSG00000228327                False  RP11-206L10.2
ENSG00000237491  ENSG00000237491                False      LINC01409
ENSG00000177757  ENSG00000177757                False         FAM87B

                feature_reference feature_biotype feature_length
index
ENSG00000241860    NCBITaxon:9606            gene           7559
ENSG00000228463    NCBITaxon:9606            gene           8224
ENSG00000228327    NCBITaxon:9606            gene           6432
ENSG00000237491    NCBITaxon:9606            gene           8413
ENSG00000177757    NCBITaxon:9606            gene           1947

adata.var.shape
(19604, 6)
```
- There are 19604 genes

- Count data is store under the X layer:

```
adata.X.A[1:10, 1:10]
array([[0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 1., 0., 0., 0., 0., 0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                            nCount_RNA  nFeature_RNA     sample  ... self_reported_ethnicity                         development_stage observation_joinid
index                                                            ...
homo17.10_w22_GTGTTGCGTTCT     22926.0          5289  homo17.10  ...                 unknown  22nd week post-fertilization human stage         8!Dr6YxLMp
homo17.09_w22_TTACGATTGTGG     16461.0          4862  homo17.09  ...                 unknown  22nd week post-fertilization human stage         jCsEC2bw18
homo17.04_w22_AAGATACCGTAC     15215.0          4155  homo17.04  ...                 unknown  22nd week post-fertilization human stage         UUAQ|A^mrh
homo17.08_w22_ATATGTCGGAAT     15033.0          4844  homo17.08  ...                 unknown  22nd week post-fertilization human stage         $v}aK)Dw^r
homo17.12_w22_GCTCAACATCCG     14490.0          4884  homo17.12  ...                 unknown  22nd week post-fertilization human stage         nGY=xT<FJK

adata.obs.shape
(118647, 28)
```
- There are 118647 cells

### Infant human neocortex
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("30b1c353-d4e1-4f34-9577-54dcc19ab4ff.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'age', 'cell_type', 'assay', 'disease', 'sex', 'tissue', 'development_stage'
    - Age: 240 day
    - Normal (no disease)
    - Assay (Drop-seq only)
    - Sex (unknown)
    - Development_stage: 7-month-old human stage
    - How many cell types? There is 1 colname for cell type: `cell_type` (3 cell types): 'glial cell', 'GABAergic neuron', 'glutamatergic neuron'
    - For tissue, there are 4 areas: 'prefrontal cortex', 'temporal lobe', 'parietal lobe', 'primary visual cortex'

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                            name  feature_is_filtered    feature_name
index
ENSG00000238009  ENSG00000238009                False    RP11-34P13.7  \
ENSG00000241860  ENSG00000241860                False   RP11-34P13.13
ENSG00000228463  ENSG00000228463                False      AP006222.1
ENSG00000230021  ENSG00000230021                False  RP11-206L10.17
ENSG00000228327  ENSG00000228327                False   RP11-206L10.2

                feature_reference feature_biotype feature_length
index
ENSG00000238009    NCBITaxon:9606            gene           3726
ENSG00000241860    NCBITaxon:9606            gene           7559
ENSG00000228463    NCBITaxon:9606            gene           8224
ENSG00000230021    NCBITaxon:9606            gene           5495
ENSG00000228327    NCBITaxon:9606            gene           6432

adata.var.shape
(24592, 6)
```
- There are 24592 genes

- Count data is store under the X layer:

```
adata.X.A[1:10, 1:10]
array([[0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 1.],
       [0., 0., 0., 1., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 1.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                          orig.ident  nCount_RNA  nFeature_RNA  lobe
index
BA10_PFC_ATGAGAGGTCAC           BA10     12518.0          4101   PFC  \
BA10_PFC_ACAATCAGGGTT           BA10     11369.0          3911   PFC
BA40_Par_CCATTATAATCC           BA40      9364.0          3611   Par
BA10_PFC_CCCGGTTCTCGA           BA10      8711.0          3088   PFC
BA41-42_Temp_TATTCGTTCGCG    BA41-42      8581.0          3651  Temp

adata.obs.shape
(51878, 28)
```
- There are 51878 cells