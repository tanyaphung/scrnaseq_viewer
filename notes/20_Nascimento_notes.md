## General information
- Paper: no associated paper
- Link: https://cellxgene.cziscience.com/collections/cae8bad0-39e9-4771-85a7-822b0e06de9f
- Raw counts downloaded on 2024-01-23
```
#Subpallial germinal zones and the developing human entorhinal cortex
wget https://datasets.cellxgene.cziscience.com/2f371e45-f65c-4c82-952f-073f823ced28.h5ad
mv 2f371e45-f65c-4c82-952f-073f823ced28.h5ad subpallial_germinal_zones.h5ad

#Interneuron maturation in the human entorhinal cortex
wget https://datasets.cellxgene.cziscience.com/70e33555-7bf9-47da-8b77-53bacb500009.h5ad
mv 70e33555-7bf9-47da-8b77-53bacb500009.h5ad interneuron_maturation.h5ad

#EC Stream Region at 14 days of age
wget https://datasets.cellxgene.cziscience.com/42beec5f-a621-4860-a27c-708c8417d5a8.h5ad
mv 42beec5f-a621-4860-a27c-708c8417d5a8.h5ad ec_stream_region.h5ad
```

## Data exploration
### Subpallial germinal zones and the developing human entorhinal cortex
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("subpallial_germinal_zones.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'author_age', 'age_group', 'region',
       'author_cell_type', 'cell_type', 'assay', 'disease',
       'organism', 'sex', 'tissue', 'self_reported_ethnicity',
       'development_stage'
    - author_age: '2y', '3y', '33d', '54d', '23GW', '14d', '13y', '27y', '79y', '51y', '50y'
    - age_group: 'Toddler (2y-3y)', 'Infant (14d-54d)', 'Fetal (23GW)', 'Teen (13y)', 'Adult (27y-79y)'
    - region: 'Postnatal EC', 'Germinal Zone', 'Migratory Stream', 'Embryonic EC'
    - assay: 10x 3' v3
    - disease: normal
    - organism: 'Homo sapiens'
    - sex: 'female', 'male'
    - tissue: 'entorhinal cortex', 'ganglionic eminence'
    - development_stage: '2-year-old human stage', '3-year-old human stage', '1-month-old human stage', '2-month-old human stage', '23rd week post-fertilization human stage', 'newborn human stage', '13-year-old human stage', '27-year-old human stage', '79-year-old human stage', '51-year-old human stage', '50-year-old human stage'
    - How many cell types?
        - `cell_type` (10 cell types): 'brain vascular cell', 'cortical interneuron', 'glutamatergic neuron', 'neural progenitor cell', 'oligodendrocyte', 'oligodendrocyte precursor cell', 'astrocyte', 'microglial cell', 'ependymal cell', 'inhibitory interneuron'
        - `author_cell_type` (11 cell types): 'Vascular Cells', 'Cortical Interneurons', 'Excitatory Neurons', 'Progenitor Cells', 'Oligodendrocytes', 'OPCs', 'Astrocytes', 'Microglia', 'Ependymal Cells', 'Non-cortical Interneurons', 'Dividing Cells'

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 vst.mean  vst.variance  vst.variance.expected  vst.variance.standardized  ...  feature_name  feature_reference feature_biotype feature_length
index                                                                                      ...                                                               
ENSG00000243485  0.000120      0.000120               0.000121                   0.991504  ...   MIR1302-2HG     NCBITaxon:9606            gene           1021
ENSG00000237613  0.000008      0.000008               0.000008                   1.000016  ...       FAM138A     NCBITaxon:9606            gene           1219
ENSG00000186092  0.000000      0.000000               0.000000                   0.000000  ...         OR4F5     NCBITaxon:9606            gene           2618
ENSG00000284733  0.000000      0.000000               0.000000                   0.000000  ...        OR4F29     NCBITaxon:9606            gene            939
ENSG00000284662  0.000000      0.000000               0.000000                   0.000000  ...        OR4F16     NCBITaxon:9606            gene            995


adata.var.shape
(21563, 10)
```
- There are 21,563 genes

- Count data is stored under the raw.X layer:

```
# the X layer shows decimals
metadata = adata.obs
metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
metadata = metadata.rename(columns = {'index': 'cell_id'})
cell_id = metadata[metadata["tissue"] == "ganglionic eminence"]["cell_id"]
adata.obs.set_index('index', inplace=True)
adata_subset = adata[cell_id]

adata_subset.X.A[1:10, 1:10]
array([[0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.9901169, 0.       , 0.       ],
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
       [0.       , 0.       , 0.       , 0.       , 0.7986996, 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ]], dtype=float32)
# the raw count is in the raw.X layer
adata_subset.raw.X.A[1:10, 1:10]
array([[0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 1., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 1., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                      nCount_RNA  nFeature_RNA sample    donor_id  ...             tissue  self_reported_ethnicity        development_stage observation_joinid
AAACCCAAGCTCATAC-1_1     35771.0          8320    H48  manaab_H48  ...  entorhinal cortex                  unknown   2-year-old human stage         LsC<fCieJZ
AAACCCAAGGTCTACT-1_1     21542.0          6828    H39  manaab_H39  ...  entorhinal cortex                  unknown   3-year-old human stage         uuL;8H~QbK
AAACCCAAGTCATGAA-1_1     22788.0          7539    H46  manaab_H46  ...  entorhinal cortex                  unknown  1-month-old human stage         %j-z_$ZLj-
AAACCCACAAATTGCC-1_1     29773.0          7590    H39  manaab_H39  ...  entorhinal cortex                  unknown   3-year-old human stage         reMKtQ(8%c
AAACCCACAACTTGGT-1_1      5963.0          2991    H46  manaab_H46  ...  entorhinal cortex                  unknown  1-month-old human stage         A`g?hei>wl


adata.obs.shape
(69174, 35)
```
- There are 69,174 cells

**IMPORTANT NOTES**
- Separate out by 2 tissues and 5 development stage

### Interneuron maturation in the human entorhinal cortex
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("interneuron_maturation.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'author_age', 'age_group', 'region', 'author_cell_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'development_stage'
    - author_age: '3y', '2y', '54d', '33d', '23GW', '14d', '13y', '27y', '79y', '51y', '50y'
    - age_group: 'Toddler (2y-3y)', 'Infant (14d-54d)', 'Fetal (23GW)', 'Teen (13y)', 'Adult (27y-79y)'
    - region: 'Postnatal EC', 'Germinal Zone', 'Migratory Stream', 'Embryonic EC'
    - assay: 10x 3' v3
    - disease: normal
    - organism: 'Homo sapiens'
    - sex: 'female', 'male'
    - tissue: 'entorhinal cortex', 'ganglionic eminence'
    - development_stage: '3-year-old human stage', '2-year-old human stage', '2-month-old human stage', '1-month-old human stage', '23rd week post-fertilization human stage', 'newborn human stage', '13-year-old human stage', '27-year-old human stage', '79-year-old human stage', '51-year-old human stage', '50-year-old human stage'
    - How many cell types?
        - `cell_type` (9 cell types): 'sst GABAergic cortical interneuron', 'forebrain neuroblast', 'caudal ganglionic eminence derived GABAergic cortical interneuron', 'pvalb GABAergic cortical interneuron', 'lamp5 GABAergic cortical interneuron', 'chandelier cell', 'vip GABAergic cortical interneuron', 'sncg GABAergic cortical interneuron', 'sst chodl GABAergic cortical interneuron'
        - `author_cell_type` (14 cell types): 'Sst', 'Immature Interneurons - MGE 1', 'Pax6', 'Pvalb', 'Lamp5 Lhx6', 'Chandelier', 'Vip', 'Lamp5', 'Immature Interneurons - CGE 1', 'Sncg', 'Immature Interneurons - CGE 2', 'Sst Chodl', 'Immature Interneurons - MGE 2', 'Immature Interneurons - Mix'

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 vst.mean  vst.variance  vst.variance.expected  vst.variance.standardized  ...  feature_name  feature_reference feature_biotype feature_length
index                                                                                      ...                                                               
ENSG00000243485  0.000195      0.000195               0.000188                   1.038307  ...   MIR1302-2HG     NCBITaxon:9606            gene           1021
ENSG00000237613  0.000000      0.000000               0.000000                   0.000000  ...       FAM138A     NCBITaxon:9606            gene           1219
ENSG00000186092  0.000000      0.000000               0.000000                   0.000000  ...         OR4F5     NCBITaxon:9606            gene           2618
ENSG00000284733  0.000000      0.000000               0.000000                   0.000000  ...        OR4F29     NCBITaxon:9606            gene            939
ENSG00000284662  0.000000      0.000000               0.000000                   0.000000  ...        OR4F16     NCBITaxon:9606            gene            995


adata.var.shape
(21563, 10)
```
- There are 21,563 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
metadata = adata.obs
metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
metadata = metadata.rename(columns = {'index': 'cell_id'})
cell_id = metadata[metadata["tissue"] == "ganglionic eminence"]["cell_id"]
adata.obs.set_index('index', inplace=True)
adata_subset = adata[cell_id]

adata_subset.X.A[1:5, 1:25]
array([[0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 2.4310644, 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       , 0.       , 0.       , 0.       ]],
      dtype=float32)
# the raw count is in the raw.X layer
adata_subset.raw.X.A[1:5, 1:25]
array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
        0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        2., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                      nCount_RNA  nFeature_RNA sample    donor_id  ...             tissue  self_reported_ethnicity        development_stage observation_joinid
index                                                              ...                                                                                       
AAACCCAAGGTCTACT-1_1     21542.0          6828    H39  manaab_H39  ...  entorhinal cortex                  unknown   3-year-old human stage         uuL;8H~QbK
AAACGAAGTGCCGTAC-1_1     22533.0          6542    H48  manaab_H48  ...  entorhinal cortex                  unknown   2-year-old human stage         ^#FFnP*Ue<
AAACGCTGTATGGAAT-1_1      1406.0           990    H39  manaab_H39  ...  entorhinal cortex                  unknown   3-year-old human stage         ZyJ;=L`Zi6
AAAGAACCATATCGGT-1_1     10174.0          4179    H29  manaab_H29  ...  entorhinal cortex                  unknown  2-month-old human stage         1Iw2QCL)CQ
AAAGAACGTACCGTGC-1_1     16126.0          5742    H29  manaab_H29  ...  entorhinal cortex                  unknown  2-month-old human stage         -<UmKQqD2m


adata.obs.shape
(20470, 32)
```
- There are 20470 cells

**IMPORTANT NOTES**
- Separate out by 2 tissues and 5 development stage

### EC Stream Region at 14 days of age
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("ec_stream_region.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'author_cell_type', 'cell_type', 'assay', 'disease',
       'organism', 'sex', 'tissue', 'development_stage'
    - author_age: '3y', '2y', '54d', '33d', '23GW', '14d', '13y', '27y', '79y', '51y', '50y'
    - age_group: 'Toddler (2y-3y)', 'Infant (14d-54d)', 'Fetal (23GW)', 'Teen (13y)', 'Adult (27y-79y)'
    - region: 'Postnatal EC', 'Germinal Zone', 'Migratory Stream', 'Embryonic EC'
    - assay: 10x 3' v3
    - disease: normal
    - organism: 'Homo sapiens'
    - sex: 'female', 'male'
    - tissue: 'entorhinal cortex'
    - development_stage: 'newborn human stage'
    - How many cell types?
        - `cell_type` (10 cell types): 'radial glial cell', 'astrocyte', 'ependymal cell', 'oligodendrocyte precursor cell', 'GABAergic interneuron', 'microglial cell', 'glutamatergic neuron', 'brain vascular cell', 'oligodendrocyte', 'leukocyte'
        - `author_cell_type` (14 cell types): 'Radial glia', 'Astrocytes', 'Ependymal cells', 'Early OPCs', 'OPCs', 'Young inhibitory neurons', 'Microglia', 'Differentiating OPCs', 'Intermediary progenitors', 'Young excitatory neurons', 'Vascular and perivascular cells', 'Oligodendrocytes', 'Immune system mix', 'Dividing Cells'

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 vst.mean  vst.variance  vst.variance.expected  vst.variance.standardized  ...  feature_name  feature_reference feature_biotype feature_length
index                                                                                      ...                                                               
ENSG00000243485  0.000263      0.000263               0.000258                   1.016613  ...   MIR1302-2HG     NCBITaxon:9606            gene           1021
ENSG00000237613  0.000000      0.000000               0.000000                   0.000000  ...       FAM138A     NCBITaxon:9606            gene           1219
ENSG00000186092  0.000000      0.000000               0.000000                   0.000000  ...         OR4F5     NCBITaxon:9606            gene           2618
ENSG00000284733  0.000000      0.000000               0.000000                   0.000000  ...        OR4F29     NCBITaxon:9606            gene            939
ENSG00000284662  0.000000      0.000000               0.000000                   0.000000  ...        OR4F16     NCBITaxon:9606            gene            995


adata.var.shape
(21563, 10)
```
- There are 21,563 genes

- Count data is store under the raw.X layer:

```
adata.X.A[1:5, 1:25]
array([[0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.93400043,
        0.        , 0.        , 0.94206774, 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.9550423 , 1.4306066 , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 1.3237723 ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ]], dtype=float32)
# the raw count is in the raw.X layer
adata.raw.X.A[1:5, 1:25]
array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 1., 2., 0.],
       [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                    nCount_RNA  nFeature_RNA    donor_id  percent.mt  ...             tissue self_reported_ethnicity    development_stage observation_joinid
AAACCCAAGGTCTTTG-1      4453.0          2412  manaab_H71    0.000000  ...  entorhinal cortex                 unknown  newborn human stage         c^Cc-<WRj8
AAACCCACACTCAGAT-1      6495.0          2922  manaab_H71    0.066151  ...  entorhinal cortex                 unknown  newborn human stage         #nZV4Fx(7w
AAACCCAGTACCAGAG-1      6688.0          3302  manaab_H71    0.421119  ...  entorhinal cortex                 unknown  newborn human stage         kFX+TCW+r-
AAACCCAGTGCGGTAA-1      3732.0          2070  manaab_H71    0.070207  ...  entorhinal cortex                 unknown  newborn human stage         ?+~)sJ5nO@
AAACCCAGTTCGGCGT-1      4665.0          2646  manaab_H71    0.000000  ...  entorhinal cortex                 unknown  newborn human stage         W^}d?n)2hh


adata.obs.shape
(7618, 26)
```
- There are 7,618 cells

**IMPORTANT NOTES**
- No separation needed