## General information
- Paper: Molecular signatures underlying neurofibrillary tangle susceptibility in Alzheimer’s disease
- Link: https://europepmc.org/article/med/35882228
- Raw counts downloaded on 2024-01-17
```
wget https://datasets.cellxgene.cziscience.com/9edbdc95-c2ef-4cdd-9b19-89bc402904fa.h5ad #this is for the excitatory neurons
wget https://datasets.cellxgene.cziscience.com/ad5f0723-101f-42f3-a31b-4319f6187ccf.h5ad #this is for the inhibitory neurons
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("9edbdc95-c2ef-4cdd-9b19-89bc402904fa.h5ad")
```
- Explore `adata.obs`: 'Cell.Types', 'cell_type', 'assay', 'disease',
       'organism', 'sex', 'tissue', 'development_stage'
    - disease: 'Alzheimer disease', 'normal' -> need to keep normal only
    - Assay ("10x 3' v2", "10x 3' v3")
    - Sex: female and male
    - Organism: Homo sapiens
    - Tissues: prefrontal cortex
    - Development stage: adult (['73-year-old human stage', '62-year-old human stage', '57-year-old human stage', '81-year-old human stage', '79-year-old human stage', '80 year-old and over human stage', '89-year-old human stage', '61-year-old human stage', '67-year-old human stage', '87-year-old human stage', '72-year-old human stage', '66-year-old human stage', '68-year-old human stage', '71-year-old human stage'])
    - How many cell types? 
        - `Cell.Types` (13 cell types): 'Ex04_RORB-GABRG1 (L4-L5)', 'Ex13_FEZF2-SEMA3D (L6b)', 'Ex03_RORB-MME (L4-L5)', 'Ex07_RORB-PCP4 (L5)', 'Ex02_CUX2-COL5A2 (L2-L4)', 'Ex01_CUX2-LAMP5 (L2-L3)', 'Ex12_FEZF2-SYT6 (L6)', 'Ex05_RORB-ADGRLA (L5)', 'Ex08_FEZF2-PCP4-ROBO3 (L5)', 'Ex06_RORB-PCP4-RPRM (L5)', 'Ex10_THEMIS (L5-6)', 'Ex11_THEMIS-NTNG2 (L6)', 'Ex09_FEZF2-ADRA1A (L5b)'
    - CONCLUSION: separate normal only
- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered feature_name feature_reference feature_biotype feature_length
feature_id
ENSG00000278915                False  AC005609.18    NCBITaxon:9606            gene            577
ENSG00000168454                False       TXNDC2    NCBITaxon:9606            gene           3392
ENSG00000139180                False       NDUFA9    NCBITaxon:9606            gene          10391
ENSG00000229177                False   AC008154.5    NCBITaxon:9606            gene            327
ENSG00000204564                False     C6orf136    NCBITaxon:9606            gene           2118

adata.var.shape
(33178, 5)
```
- There are 33,178 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata_subset.X.A[1:20,1:20]
array([[0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.9595668 , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.9595668 ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 1.0698613 ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.8361783 , 0.        , 0.8361783 , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.78652054,
        1.2212683 , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.78652054],
       [0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.8541368 , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        ],
       [0.        , 0.        , 0.        , 0.        , 1.3116333 ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 0.        , 0.        ,
        0.        , 0.        , 0.        , 1.3116333 ],
# the raw count is in the raw.X layer
adata_subset.raw.X.A[1:20,1:20]
array([[0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 1.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 1.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 2., 0., 0., 0., 0., 0.,
        0., 0., 1.],
       [0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 1.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
        0., 0., 1.],
       [0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
        0., 0., 2.],
       [0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0.],
       [0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.,
        0., 0., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                          nCount_RNA  nFeature_RNA  ...        development_stage observation_joinid
index                                               ...
C0001_AAACCTGAGGAGTTGC-1      2079.0          1386  ...  73-year-old human stage         J2Ev2aON(S
C0001_AAACCTGTCAATCTCT-1      7545.0          3415  ...  73-year-old human stage         4t!~r-e2M@
C0001_AAACCTGTCGTCACGG-1      3440.0          1981  ...  73-year-old human stage         11Z{9>%+ln
C0001_AAACGGGAGGAGTTTA-1      4929.0          2752  ...  73-year-old human stage         w0TZ1oDU9v
C0001_AAACGGGAGTGCCATT-1      5467.0          2823  ...  73-year-old human stage         ecYuVdMWa6

[5 rows x 36 columns]

adata.obs.shape
(96129, 36)
```
- For excitatory neurons, there are 96,129 cells total and 46,070 are normal
- For inhibitory neurons, there are 23,197 cells total and 11,464 are normal
- Total: 57,534 cells
- CONCLUSIONS:
    - 2 files: one for excitatory neuron and one for inhibitory neuron
    - Keep normal only
- There are 2 files downloaded
```
9edbdc95-c2ef-4cdd-9b19-89bc402904fa.h5ad
ad5f0723-101f-42f3-a31b-4319f6187ccf.h5ad
```

But the donor id are the same so I'm going to combine into 1 file after the step to compute mean and specifity for magma gene property analyses. 
```
import anndata
adata = anndata.read("9edbdc95-c2ef-4cdd-9b19-89bc402904fa.h5ad")
for i in adata.obs['donor_id'].unique():
...     print(i)
...
Subject6
Subject8
Subject4
Subject5
Subject2
Subject3
Subject1
Subject7
CTRL-1
CTRL-2
CTRL-3
CTRL-4
CTRL-5
CTRL-6
CTRL-7
CTRL-8

adata = anndata.read("ad5f0723-101f-42f3-a31b-4319f6187ccf.h5ad")
for i in adata.obs['donor_id'].unique():
...     print(i)
...
Subject6
Subject8
Subject4
Subject5
Subject2
Subject3
Subject1
Subject7
CTRL-1
CTRL-2
CTRL-3
CTRL-4
CTRL-5
CTRL-6
CTRL-7
CTRL-8
```