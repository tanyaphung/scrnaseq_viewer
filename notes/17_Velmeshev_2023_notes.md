## General information
- Raw counts downloaded on 2023-12-19
```
wget https://datasets.cellxgene.cziscience.com/e9e38c8d-fd21-475d-acf6-f8994b003d70.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("e9e38c8d-fd21-475d-acf6-f8994b003d70.h5ad")
```

- What are the developmental stages: 

```
for i in list(adata.obs["development_stage"].unique()):
...     print(i)
...
12-year-old human stage
14-year-old human stage
22-year-old human stage
21-year-old human stage
19-year-old human stage
13-year-old human stage
15-year-old human stage
8-year-old human stage
6-year-old human stage
4-year-old human stage
39-year-old human stage
54-year-old human stage
34-year-old human stage
44-year-old human stage
2-year-old human stage
under-1-year-old human stage
ninth LMP month human stage
2-month-old human stage
3-year-old human stage
1-month-old human stage
1-year-old human stage
fifth LMP month human stage
eighth LMP month human stage
sixth LMP month human stage
seventh LMP month human stage
5-year-old human stage
3-month-old human stage
16-year-old human stage
fourth LMP month human stage
45-year-old human stage
53-year-old human stage
28-year-old human stage
newborn human stage
4-month-old human stage
5-month-old human stage
10-year-old human stage
20-year-old human stage
40-year-old human stage
17-year-old human stage
25-year-old human stage
```
- I will group them into 8 groups as in Figure 1 from the paper (https://www.science.org/doi/10.1126/science.adf0834#sec-1):
  
| developmental stage | definition | group | 
| ------------------- | ---------- | ----- | 
| 12-year-old human stage | 10-20 years | group_7 | 
| 14-year-old human stage | 10-20 years | group_7 | 
| 22-year-old human stage | Adult | group_8 |
| 21-year-old human stage | Adult | group_8 |
| 19-year-old human stage | 10-20 years | group_7 |
| 13-year-old human stage | 10-20 years | group_7 |
| 15-year-old human stage | 10-20 years | group_7 |
| 8-year-old human stage | 4-10 years | group_6 | 
| 6-year-old human stage | 4-10 years | group_6 | 
| 4-year-old human stage | 4-10 years | group_6 | 
| 39-year-old human stage | Adult | group_8 |
| 54-year-old human stage | Adult | group_8 |
| 34-year-old human stage | Adult | group_8 |
| 44-year-old human stage | Adult | group_8 |
| 2-year-old human stage | 2-4 years | group_5 |
| under-1-year-old human stage | 0-1 years | group_3 |
| ninth LMP month human stage | 3rd trimester | group_2 |
| 2-month-old human stage | 0-1 years | group_3 |
| 3-year-old human stage | 2-4 years | group_5 |
| 1-month-old human stage | 0-1 years | group_3 |
| 1-year-old human stage | 1-2 years | group_4 |
| fifth LMP month human stage | 2nd trimester | group_1 |
| eighth LMP month human stage | 3rd trimester | group_2 |
| sixth LMP month human stage | 2nd trimester | group_1 |
| seventh LMP month human stage | 3rd trimester | group_2 |
| 5-year-old human stage | 4-10 years | group_6 | 
| 3-month-old human stage | 0-1 years | group_3 |
| 16-year-old human stage | 10-20 years | group_7 |
| fourth LMP month human stage | 2nd trimester | group_1 |
| 45-year-old human stage | Adult | group_8 |
| 53-year-old human stage | Adult | group_8 |
| 28-year-old human stage | Adult | group_8 |
| newborn human stage | 0-1 years | group_3 |
| 4-month-old human stage | 0-1 years | group_3 |
| 5-month-old human stage | 0-1 years | group_3 |
| 10-year-old human stage | 10-20 years | group_7 |
| 20-year-old human stage | Adult | group_8 |
| 40-year-old human stage | Adult | group_8 |
| 17-year-old human stage | 10-20 years | group_7 |
| 25-year-old human stage | Adult | group_8 |


- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered feature_name feature_reference
ENSG00000234661                False     CHL1-AS1    NCBITaxon:9606  \
ENSG00000154764                False        WNT7A    NCBITaxon:9606
ENSG00000134121                False         CHL1    NCBITaxon:9606
ENSG00000227110                False    LMCD1-AS1    NCBITaxon:9606
ENSG00000144635                False     DYNC1LI1    NCBITaxon:9606

                feature_biotype feature_length
ENSG00000234661            gene           1149
ENSG00000154764            gene           4529
ENSG00000134121            gene          11694
ENSG00000227110            gene          11257
ENSG00000144635            gene           6113

adata.var.shape
adata.var.shape
(17663, 5)
```
- There are 17663 genes

- Count data is store under the raw layer: 

```
adata.raw.X.A
array([[0., 0., 0., ..., 0., 1., 0.],
       [0., 0., 0., ..., 0., 0., 0.],
       [0., 0., 2., ..., 0., 0., 0.],
       ...,
       [0., 0., 2., ..., 0., 0., 0.],
       [0., 0., 1., ..., 0., 3., 0.],
       [0., 0., 5., ..., 0., 2., 1.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                organism_ontology_term_id
U01_AAACCTGAGAAACCAT-1_5387_BA9            NCBITaxon:9606  \
U01_AAACCTGAGACTCGGA-1_5387_BA9            NCBITaxon:9606
U01_AAACCTGCAGGGATTG-1_5387_BA9            NCBITaxon:9606
U01_AAACCTGCAGGTCGTC-1_5387_BA9            NCBITaxon:9606
U01_AAACCTGGTACGCTGC-1_5387_BA9            NCBITaxon:9606

                                tissue_ontology_term_id
U01_AAACCTGAGAAACCAT-1_5387_BA9          UBERON:0001870  \
U01_AAACCTGAGACTCGGA-1_5387_BA9          UBERON:0001870
U01_AAACCTGCAGGGATTG-1_5387_BA9          UBERON:0001870
U01_AAACCTGCAGGTCGTC-1_5387_BA9          UBERON:0001870
U01_AAACCTGGTACGCTGC-1_5387_BA9          UBERON:0001870

adata.obs.shape
(709372, 25)
```
- There are 709372 cells

- Cell type annotation
    - Just 1 level, store in `cell_type`
    
    ```
    adata.obs["cell_type"].unique()
    ['oligodendrocyte precursor cell', 'native cell', 'neural cell', 'oligodendrocyte', 'microglial cell', 'astrocyte']
    Categories (6, object): ['native cell', 'astrocyte', 'oligodendrocyte', 'microglial cell',
                            'neural cell', 'oligodendrocyte precursor cell']
    ``` 
