## General information
- Paper: Molecular and cellular evolution of the primate dorsolateral prefrontal cortex
- Link: https://www.science.org/doi/10.1126/science.abo7257
- Raw counts downloaded on 2024-01-16

```
wget https://datasets.cellxgene.cziscience.com/2e0b82ee-5320-437e-b06e-3d2150e70395.h5ad
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
- Explore `adata.obs`: 'subtype', 'subclass', 'class', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'observation_joinid'
    - Normal (no disease)
    - Assay ('10x 3' v3')
    - Sex: female and male
    - Organism: 'Homo sapiens', 'Pan troglodytes', 'Macaca mulatta', 'Callithrix jacchus'. Need to subset for human only
    - Tissues: dorsolateral prefrontal cortex
    - Development stage: adult ('64-year-old human stage', '36-year-old human stage', '19-year-old human stage', '50-year-old human stage', 'prime adult stage')
    - How many cell types? There are 3 colnames for cell type: 
        - `class` (4 cell types)
        - `subclass` (29 cell types)
        - `subtype` (114 cell types)
    - CONCLUSION: separate human only
- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                         name  gene_symbols  feature_is_filtered
ENSG00000272438  RP11-54O7.16  RP11-54O7.16                False  \
ENSG00000230699   RP11-54O7.1   RP11-54O7.1                False
ENSG00000241180   RP11-54O7.2   RP11-54O7.2                False
ENSG00000223764   RP11-54O7.3   RP11-54O7.3                False
ENSG00000187634        SAMD11        SAMD11                False

                 feature_name feature_reference feature_biotype feature_length
ENSG00000272438  RP11-54O7.16    NCBITaxon:9606            gene            351
ENSG00000230699   RP11-54O7.1    NCBITaxon:9606            gene           3043
ENSG00000241180   RP11-54O7.2    NCBITaxon:9606            gene            443
ENSG00000223764     LINC02593    NCBITaxon:9606            gene           4152
ENSG00000187634        SAMD11    NCBITaxon:9606            gene           4172

adata.var.shape
(27983, 9)
```
- There are 27,983 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata.X.A[1:25, 1:25]
[0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  1.2085186  0.         0.         0.         0.         0.
  0.         0.         0.         1.2085186  0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  1.9140283  0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.6314628  0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.        ]
# the raw count is in the raw.X layer
adata.raw.X.A[1:30, 1:30]
[[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.
  0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
  0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
:
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                           nGene     nUMI  percent.mt  mapped_reads
HSB106_1_AAACCCAAGAGGCGTT    884   1260.0    0.030479        3281.0  \
HSB106_1_AAACCCAAGGTAGTCG   1093   1651.0    0.415858        3607.0
HSB106_1_AAACCCAAGGTTATAG   1062   1832.0    0.502661        3382.0
HSB106_1_AAACCCAAGTCGCCCA   2877   7085.0    0.222627       14823.0
HSB106_1_AAACCCACAAAGGAGA   4917  13925.0    0.632122       36860.0

                                            subtype   subclass class
HSB106_1_AAACCCAAGAGGCGTT          Oligo MOG OPALIN      Oligo  Glia  \
HSB106_1_AAACCCAAGGTAGTCG          Oligo MOG OPALIN      Oligo  Glia
HSB106_1_AAACCCAAGGTTATAG          Oligo MOG OPALIN      Oligo  Glia
HSB106_1_AAACCCAAGTCGCCCA  L3-5 RORB PCBP3 IL1RAPL2  L3-5 IT-1   ExN
HSB106_1_AAACCCACAAAGGAGA   L2-3 CUX2 ACVR1C THSD7A    L2-3 IT   ExN

adata.obs.shape
(610719, 34)
```
- There are 610719 cells total and 172120 are human