## General information
- Raw counts downloaded on 2023-12-21
```
wget https://datasets.cellxgene.cziscience.com/6b141abe-83b3-439d-a66d-a3b5495c52c7.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("6b141abe-83b3-439d-a66d-a3b5495c52c7.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: cell_type, assay, disease, tissue, author_stage
    - For developmental stage, there are 18 unique values (column author_stage). However, I am going to group these to form 11 groups:
        - CS22 (CS22)
        - 11 wpc (11 wpc)
        - newborn (newborn)
        - infant (9 months) (infant)
        - infant (6 months) (infant)
        - newborn (5 days) (newborn)
        - toddler (3.5 years) (toddler)
        - 9 wpc (9 wpc)
        - 17 wpc (17 wpc)
        - adult (46 years) (adult)
        - CS18 (CS18)
        - infant (7 months) (infant)
        - adult (42 years) (adult)
        - adult (52 years) (adult)
        - adult (44 years) (adult)
        - CS19 (CS19)
        - toddler (2.8 years) (toddler)
        - 20 wpc (20 wpc)
         
        | developmental stage | definition | group |
        | ------------------- | ---------- | ----- |
        | CS18 | CS18 | group_1 |
        | CS19 | CS19 | group_2 |
        | CS22 | CS22 | group_3 |
        | 9 wpc | 9 wpc | group_4 |
        | 11 wpc | 11 wpc | group_5 |
        | 17 wpc | 17 wpc | group_6 |
        | 20 wpc | 20 wpc | group_7 |
        | newborn | newborn | group_8 |
        | infant | infant | group_9 |
        | toddler | toddler | group_10 |
        | adult | adult | group_11 |


    - How many cell types? The authors annotated 25 cell types
        - native cell
        - cerebellar granule cell
        - glutamatergic neuron
        - macroglial cell
        - Purkinje cell
        - neuroblast (sensu Vertebrata)
        - GABAergic neuron
        - cerebellar granule cell precursor
        - interneuron
        - unipolar brush cell
        - microglial cell
        - progenitor cell
        - glioblast
        - noradrenergic cell
        - central nervous system macrophage
        - brain vascular cell
        - erythroid lineage cell
        - leukocyte
        - oligodendrocyte precursor cell
        - T cell
        - meningeal macrophage
        - Bergmann glial cell
        - immature astrocyte
        - oligodendrocyte
        - differentiation-committed oligodendrocyte precursor

    - For assay, there seems to be 2 versions of 10x: 10x 3' v2, 10x 3' v3. For now I won't separate them but keep in mind there might be some differences. See Extended Data Fig. 1. 
    - For tissue, there are 5 areas but for now I won't separate them and just refer to them as cerebellum: cerebellum, cerebellar cortex, dentate nucleus, cerebellum vermis lobule, hemisphere part of cerebellar posterior lobe
**How to filter/sepeate the data?**
1. "Normal" only
2. 


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

- Count data is store under the X layer:

```
print(adata.X.A[1:10,1:10])
[[0. 0. 0. 0. 0. 0. 0. 2. 0.]
 [2. 0. 0. 0. 0. 1. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 9. 0.]
 [1. 0. 0. 0. 0. 1. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [2. 0. 0. 0. 0. 0. 0. 1. 0.]
 [1. 0. 0. 0. 0. 0. 0. 1. 0.]
 [2. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 1. 0. 0. 0.]]
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
