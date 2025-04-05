## General information
- Paper: Cryptic exon detection and transcriptomic changes revealed in single-nuclei RNA sequencing of C9ORF72 patients spanning the ALS-FTD spectrum
- Link: https://link.springer.com/article/10.1007/s00401-023-02599-5
- Raw counts downloaded on 2024-01-22
```
wget https://datasets.cellxgene.cziscience.com/cc6419c4-86ec-484a-8ff6-57c4f6675bed.h5ad

wget https://datasets.cellxgene.cziscience.com/b816a47c-1736-46c6-b54e-0217d2a10d50.h5ad

wget https://datasets.cellxgene.cziscience.com/e8acf5ae-600c-4794-a214-f3b612c82ebe.h5ad

wget https://datasets.cellxgene.cziscience.com/4ff9a4d5-7492-4c62-9c54-a5929d72f508.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("4ff9a4d5-7492-4c62-9c54-a5929d72f508.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'cluster_layers', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'development_stage'
    - assay: 10x 3' v3
    - disease: 'amyotrophic lateral sclerosis', 'amyotrophic lateral sclerosis 26 with or without frontotemporal dementia', 'normal' --> SUBSET for normal only
    - organism: 'Homo sapiens'
    - sex: 'male', 'female'
    - tissue: 'frontal cortex'
    - development_stage: '66-year-old human stage', '70-year-old human stage', '73-year-old human stage', '47-year-old human stage', '72-year-old human stage', '59-year-old human stage', '54-year-old human stage', '68-year-old human stage', '61-year-old human stage', '43-year-old human stage', '71-year-old human stage', '56-year-old human stage', '58-year-old human stage', '51-year-old human stage', '50-year-old human stage'
    - How many cell types?
        - `cell_type` (6 cell types): 'astrocyte', 'neuron', 'oligodendrocyte', 'oligodendrocyte precursor cell', 'microglial cell', 'endothelial cell'
        - `cluster_layers` (17 cell types): 'Astrocytes', 'L3-L5 Intratelencephalic Type 1', 'Oligodendrocytes', 'L2-L3 Intratelencephalic', 'OPCs', 'Somatostatin Interneurons', 'L5-L6 Near Projecting', 'Microglia', 'L3-L5 Intratelencephalic Type 2', 'L6 Intratelencephalic - Type 1', 'Endothelial', 'VIP Interneurons', 'SV2C LAMP5 Interneurons', 'L6 Corticothalamic / L6B', 'Parvalbumin interneurons', 'L6 Intratelencephalic - Type 2', 'L5 Extratelencephalic'

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                    gene_ids    feature_types     mt   ribo     hb  ...  feature_is_filtered   feature_name  feature_reference  feature_biotype  feature_length
feature_id                                                          ...                                                                                      
ENSG00000243485  MIR1302-2HG  Gene Expression  False  False  False  ...                False    MIR1302-2HG     NCBITaxon:9606             gene            1021
ENSG00000186092        OR4F5  Gene Expression  False  False  False  ...                False          OR4F5     NCBITaxon:9606             gene            2618
ENSG00000238009   AL627309.1  Gene Expression  False  False  False  ...                False   RP11-34P13.7     NCBITaxon:9606             gene            3726
ENSG00000239945   AL627309.3  Gene Expression  False  False  False  ...                False   RP11-34P13.8     NCBITaxon:9606             gene            1319
ENSG00000239906   AL627309.2  Gene Expression  False  False  False  ...                False  RP11-34P13.14     NCBITaxon:9606             gene             323


adata.var.shape
(34568, 19)
```
- There are 34,568 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
metadata = adata.obs
metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
metadata = metadata.rename(columns = {'index': 'cell_id'})
cell_id = metadata[metadata["disease"] == "normal"]["cell_id"]
adata.obs.set_index('index', inplace=True)
adata_subset = adata[cell_id]

adata_subset.X.A[1:10, 1:10]
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
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 1.2504445, 0.       ]], dtype=float32)
# the raw count is in the raw.X layer
adata_subset.raw.X.A[1:10, 1:10]
array([[0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 1., 0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
       orig.ident  nCount_RNA  nFeature_RNA  percent.mito   S.Score
V3_1          01k       617.0           482      0.012966 -0.012428  \
V224_1        01k      1733.0          1147      0.010964 -0.034455
V364_1        01k      1475.0           974      0.027797 -0.026429
V372_1        01k       966.0           718      0.010352  0.060085
V434_1        01k      1266.0           977      0.030806 -0.029667

adata.obs.shape
(69174, 35)
```

**IMPORTANT NOTES**
- https://cellxgene.cziscience.com/collections/aee9c366-f2fb-470b-8937-577d5d87d3fc
- There are 4 files: 
    - 341_Gittings_FrontalCortex_Human_2023_Part1 & 342_Gittings_FrontalCortex_Human_2023_Part2
        - In this one, the index which is the cell id is `barcode`
        - In this one, there cell types columns are: `cell_type` and `cluster_layers`
    - 343_Gittings_OccipitalCortex_Human_2023_Part1 & 344_Gittings_OccipitalCortex_Human_2023_Part2
        - the index is not label so it's just `index`
        - cell type columns are `cell_type` and `cluster`