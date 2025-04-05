## General information
- Paper: Comparative cellular analysis of motor cortex in human, marmoset and mouse
- Link: https://www.nature.com/articles/s41586-021-03465-8#Abs1
- Note: the way that this data is available for download is a bit strange. There are 4 files available for download:
    - Evolution of cellular diversity in primary motor cortex of human, marmoset monkey, and mouse 3-species integration inhibitory neurons
    - Evolution of cellular diversity in primary motor cortex of human, marmoset monkey, and mouse 4-species integration excitory neurons
    - Evolution of cellular diversity in primary motor cortex of human, marmoset monkey, and mouse 3-species integration excitory neurons
    - Evolution of cellular diversity in primary motor cortex of human, marmoset monkey, and mouse 3-species integration non-nuerons
    - It seemns that there are 3 files, each for inhibitory neurons, excitatory neurons, and non-neurons. These 3 files were integrated together with 3 species. In addition, there is 1 file for excitatory neurons that were integrated together with 4 species.
- Raw counts downloaded on 2024-01-08
```
# Inhibitory neurons 3 species integration
wget https://datasets.cellxgene.cziscience.com/62ab3ab7-b0a7-409e-9ed6-4a50d60b964a.h5ad

# Excitatory neurons 4 species integration
wget https://datasets.cellxgene.cziscience.com/d0e7effe-2311-4c69-93d5-462ffb6bc663.h5ad

# Excitarory neurons 3 species integration
wget https://datasets.cellxgene.cziscience.com/f8d20729-9133-4113-b49a-3f996e5f8a78.h5ad

# Non-neurons 3 species integration
wget https://datasets.cellxgene.cziscience.com/c0c1469a-1696-450d-949c-7c3a13299ed3.h5ad
```

## Data exploration
### inhibitory_neurons_3speciesintegration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("62ab3ab7-b0a7-409e-9ed6-4a50d60b964a.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'donor_id', 'BICCN_cluster_label', 'BICCN_subclass_label', 'subclass_id', 'sample_id', 'cell_type', 'assay', 'disease', 'organism', 'sex',
       'tissue', 'development_stage'
    - Normal (no disease)
    - Assay (10x 3' v3)
    - Sex: female and male
    - Organism: 'Callithrix jacchus', 'Homo sapiens', 'Mus musculus'. Need to subset for human only
    - Tissue: primary motor cortex
    - Development stage: adult ('50-year-old human stage', '60-year-old human stage', 'early adult stage', 'post-juvenile adult stage')
    - How many cell types? There are 3 colnames for cell type: 
        - `cell_type` (1 cell type): GABAergic neuron
        - `BICCN_subclass_label` (7 cell types): 'Sst', 'Pvalb', 'Vip', 'Lamp5', 'Sncg', 'Sst Chodl', 'Meis2'
        - `BICCN_cluster_label` (183 cell types)

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered feature_name feature_reference feature_biotype feature_length
feature_id
ENSG00000140463                False         BBS4    NCBITaxon:9606            gene           4148
ENSG00000125533                False      BHLHE23    NCBITaxon:9606            gene           1109
ENSG00000182667                False          NTM    NCBITaxon:9606            gene           6956
ENSG00000101440                False         ASIP    NCBITaxon:9606            gene            845
ENSG00000100568                False        VTI1B    NCBITaxon:9606            gene           5768

adata.var.shape
(14736, 5)
```
- There are 14,736 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata.X.A[1:10, 1:10]
array([[0.       , 2.3978953, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.6931472],
       [0.       , 1.7917595, 0.       , 0.6931472, 0.       , 0.       ,
        0.6931472, 0.       , 0.       ],
       [0.       , 2.6390574, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 1.609438 , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.6931472, 0.       ],
       [0.       , 3.0445225, 0.       , 0.6931472, 0.       , 0.       ,
        0.6931472, 0.       , 0.       ],
       [0.       , 3.912023 , 0.       , 0.       , 0.       , 0.       ,
        0.6931472, 0.       , 0.       ],
       [0.       , 3.6635616, 0.       , 0.       , 0.       , 0.       ,
        1.3862944, 0.       , 0.       ],
       [0.       , 2.8903718, 0.       , 0.       , 0.       , 0.6931472,
        0.       , 0.       , 0.       ],
       [0.       , 2.7725887, 0.       , 0.       , 0.       , 0.       ,
        0.6931472, 0.       , 0.       ]], dtype=float32)
# the raw count is in the raw.X layer
adata.raw.X.A[1:10, 1:10]
array([[  0.,   8.,   0.,   0.,   0.,   0.,   0.,   0.,   1.],
       [  0.,   4.,   0.,   1.,   0.,   0.,   1.,   0.,   0.],
       [  0.,  11.,   0.,   0.,   0.,   0.,   0.,   0.,   0.],
       [  0.,   3.,   0.,   0.,   0.,   0.,   0.,   1.,   0.],
       [  0.,  22.,   0.,   1.,   0.,   0.,   1.,   0.,   0.],
       [  0.,  62.,   0.,   0.,   0.,   0.,   1.,   0.,   0.],
       [  0., 122.,   0.,   0.,   0.,   0.,   8.,   0.,   0.],
       [  0.,  23.,   0.,   0.,   0.,   1.,   0.,   0.,   0.],
       [  0.,  11.,   0.,   0.,   0.,   0.,   1.,   0.,   0.]],
      dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                              nCount_RNA  nFeature_RNA  ...        development_stage observation_joinid
human_AAACCCAAGGATTTCC-21L8TX_180927_001_A01     12053.0          3994  ...  60-year-old human stage         qHMCh(28^i
human_AAAGAACAGGTGAGCT-21L8TX_180927_001_A01      9523.0          3526  ...  60-year-old human stage         KMtgvNUQQ|
human_AAAGAACGTACAGCGA-21L8TX_180927_001_A01      9177.0          3727  ...  60-year-old human stage         jH)C|l@9n|
human_AAATGGAGTGAGTTTC-21L8TX_180927_001_A01     10083.0          3975  ...  60-year-old human stage         A<arUTit=w
human_AACAAAGCAACAAGTA-21L8TX_180927_001_A01      9956.0          3424  ...  60-year-old human stage         rT=S;nqQ#l

adata.obs.shape
(29486, 36)
```
- There are 29,486 cells

### excitatory_neurons_4speciesintegration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("d0e7effe-2311-4c69-93d5-462ffb6bc663.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'sample_id', 'BICCN_subclass_label', 'BICCN_cluster_label', 'donor_id','cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'development_stage'
    - Normal (no disease)
    - Assay (10x 3' v3)
    - Sex: female, male, and unknown
    - Organism: 'Callithrix jacchus', 'Macaca mulatta', 'Homo sapiens', 'Mus musculus'. Need to subset for human only
    - Tissue: primary motor cortex
    - Development stage: adult ('50-year-old human stage', '60-year-old human stage', 'early adult stage', 'post-juvenile adult stage')
    - How many cell types? There are 3 colnames for cell type: 
        - `cell_type` (1 cell type): glutamatergic neuron
        - `BICCN_subclass_label` (8 cell types): 'L2/3 IT', 'L5 ET', 'L5 IT', 'L5/6 NP', 'L6 CT', 'L6 IT', 'L6 IT Car3', 'L6b'
        - `BICCN_cluster_label` (112 cell types)

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered feature_name feature_reference feature_biotype feature_length
feature_id
ENSG00000140463                False         BBS4    NCBITaxon:9606            gene           4148
ENSG00000125533                False      BHLHE23    NCBITaxon:9606            gene           1109
ENSG00000182667                False          NTM    NCBITaxon:9606            gene           6956
ENSG00000101440                False         ASIP    NCBITaxon:9606            gene            845
ENSG00000100568                False        VTI1B    NCBITaxon:9606            gene           5768

adata.var.shape
(14617, 5)
```
- There are 14,617 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata.X.A[1:10, 1:10]
array([[0.       , 0.       , 0.       , 0.6931472, 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 1.0986123, 0.       , 0.       , 0.       , 0.       ,
        1.3862944, 0.       , 0.       ],
       [0.       , 3.218876 , 0.       , 0.       , 0.       , 0.       ,
        0.6931472, 0.       , 0.       ],
       [0.       , 1.7917595, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 2.6390574, 0.       , 0.6931472, 0.       , 0.       ,
        0.6931472, 0.       , 0.       ],
       [0.       , 4.26268  , 0.       , 0.       , 0.       , 0.       ,
        1.0986123, 0.       , 0.       ],
       [0.       , 2.8332133, 0.       , 0.       , 0.       , 0.6931472,
        0.6931472, 0.       , 0.       ],
       [0.       , 2.0794415, 0.       , 0.       , 0.       , 0.       ,
        1.0986123, 0.       , 0.       ],
       [0.       , 3.713572 , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ]], dtype=float32)
# the raw count is in the raw.X layer
adata.raw.X.A[1:10, 1:10]
array([[  0.,   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.],
       [  0.,   0.,   0.,   0.,   0.,   0.,   2.,   0.,   0.],
       [  0.,  29.,   0.,   0.,   0.,   0.,   1.,   0.,   0.],
       [  0.,   4.,   0.,   0.,   0.,   0.,   0.,   0.,   0.],
       [  0.,  12.,   0.,   1.,   0.,   0.,   1.,   0.,   0.],
       [  0., 108.,   0.,   0.,   0.,   0.,   3.,   0.,   0.],
       [  0.,  19.,   0.,   0.,   0.,   1.,   1.,   0.,   0.],
       [  0.,  18.,   0.,   1.,   0.,   0.,   4.,   0.,   0.],
       [  0.,  21.,   0.,   0.,   0.,   0.,   0.,   0.,   0.]],
      dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                              nCount_RNA  nFeature_RNA  ...        development_stage observation_joinid
human_AAACCCAAGTATGGCG-21L8TX_180927_001_A01      9612.0          3563  ...  60-year-old human stage         TgP04NBY-G
human_AAACGAAGTATGAAGT-21L8TX_180927_001_A01     25235.0          5652  ...  60-year-old human stage         $A^KJkuQf=
human_AAACGCTAGTGCACCC-21L8TX_180927_001_A01      9890.0          3716  ...  60-year-old human stage         aKC#T949qp
human_AAAGGATAGAATTCAG-21L8TX_180927_001_A01     25754.0          5841  ...  60-year-old human stage         %=6UF@}!a2
human_AAAGGTAAGGGCGAGA-21L8TX_180927_001_A01     16644.0          4573  ...  60-year-old human stage         9PLyC;$we+

adata.obs.shape
(29050, 35)
```
- There are 29,050 cells

### excitatory_neurons_3speciesintegration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("f8d20729-9133-4113-b49a-3f996e5f8a78.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'donor_id', 'BICCN_cluster_label', 'BICCN_subclass_label', 'subclass_id', 'sample_id', 'cell_type', 'assay', 'disease', 'organism', 'sex',
       'tissue', 'development_stage'
    - Normal (no disease)
    - Assay (10x 3' v3)
    - Sex: female, male
    - Organism: 'Callithrix jacchus', 'Homo sapiens', 'Mus musculus'. Need to subset for human only
    - Tissue: primary motor cortex
    - Development stage: adult ('50-year-old human stage', '60-year-old human stage', 'early adult stage', 'post-juvenile adult stage')
    - How many cell types? There is 2 colname for cell type: 
        - `cell_type` (1 cell type): glutamatergic neuron
        - `BICCN_subclass_label` (8 cell types): 'L2/3 IT', 'L5 ET', 'L5 IT', 'L5/6 NP', 'L6 CT', 'L6 IT', 'L6 IT Car3', 'L6b'
        - `BICCN_cluster_label` (99 cell types)

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered feature_name feature_reference feature_biotype feature_length
feature_id
ENSG00000140463                False         BBS4    NCBITaxon:9606            gene           4148
ENSG00000125533                False      BHLHE23    NCBITaxon:9606            gene           1109
ENSG00000182667                False          NTM    NCBITaxon:9606            gene           6956
ENSG00000101440                False         ASIP    NCBITaxon:9606            gene            845
ENSG00000100568                False        VTI1B    NCBITaxon:9606            gene           5768

adata.var.shape
(14763, 5)
```
- There are 14,763 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata.X.A[1:10, 1:10]
array([[0.       , 0.       , 0.       , 0.6931472, 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 1.0986123, 0.       , 0.       , 0.       , 0.       ,
        1.3862944, 0.       , 0.       ],
       [0.       , 3.218876 , 0.       , 0.       , 0.       , 0.       ,
        0.6931472, 0.       , 0.       ],
       [0.       , 1.7917595, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 2.6390574, 0.       , 0.6931472, 0.       , 0.       ,
        0.6931472, 0.       , 0.       ],
       [0.       , 4.26268  , 0.       , 0.       , 0.       , 0.       ,
        1.0986123, 0.       , 0.       ],
       [0.       , 2.8332133, 0.       , 0.       , 0.       , 0.6931472,
        0.6931472, 0.       , 0.       ],
       [0.       , 2.0794415, 0.       , 0.       , 0.       , 0.       ,
        1.0986123, 0.       , 0.       ],
       [0.       , 3.713572 , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ]], dtype=float32)
# the raw count is in the raw.X layer
adata.raw.X.A[1:10, 1:10]
array([[  0.,   0.,   0.,   1.,   0.,   0.,   0.,   0.,   0.],
       [  0.,   0.,   0.,   0.,   0.,   0.,   2.,   0.,   0.],
       [  0.,  29.,   0.,   0.,   0.,   0.,   1.,   0.,   0.],
       [  0.,   4.,   0.,   0.,   0.,   0.,   0.,   0.,   0.],
       [  0.,  12.,   0.,   1.,   0.,   0.,   1.,   0.,   0.],
       [  0., 108.,   0.,   0.,   0.,   0.,   3.,   0.,   0.],
       [  0.,  19.,   0.,   0.,   0.,   1.,   1.,   0.,   0.],
       [  0.,  18.,   0.,   1.,   0.,   0.,   4.,   0.,   0.],
       [  0.,  21.,   0.,   0.,   0.,   0.,   0.,   0.,   0.]],
      dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                              nCount_RNA  nFeature_RNA  ...        development_stage observation_joinid
human_AAACCCAAGTATGGCG-21L8TX_180927_001_A01      9670.0          3579  ...  60-year-old human stage         TgP04NBY-G
human_AAACGAAGTATGAAGT-21L8TX_180927_001_A01     25409.0          5697  ...  60-year-old human stage         $A^KJkuQf=
human_AAACGCTAGTGCACCC-21L8TX_180927_001_A01      9948.0          3736  ...  60-year-old human stage         aKC#T949qp
human_AAAGGATAGAATTCAG-21L8TX_180927_001_A01     25874.0          5876  ...  60-year-old human stage         %=6UF@}!a2
human_AAAGGTAAGGGCGAGA-21L8TX_180927_001_A01     16732.0          4604  ...  60-year-old human stage         9PLyC;$we+

adata.obs.shape
(24213, 36)
```
- There are 24,213 cells

### non_neurons_3speciesintegration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("c0c1469a-1696-450d-949c-7c3a13299ed3.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'donor_id', 'BICCN_cluster_label', 'BICCN_subclass_label', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue',  'development_stage'
    - Normal (no disease)
    - Assay (10x 3' v3)
    - Sex: female, male
    - Organism: 'Callithrix jacchus', 'Homo sapiens', 'Mus musculus'. Need to subset for human only
    - Tissue: primary motor cortex
    - Development stage: adult ('50-year-old human stage', '60-year-old human stage', 'early adult stage', 'post-juvenile adult stage')
    - How many cell types? There is 2 colname for cell type: 
        - `cell_type` (8 cell type): 'oligodendrocyte', 'microglial cell', 'oligodendrocyte precursor cell', 'astrocyte', 'endothelial cell', 'leptomeningeal cell', 'pericyte', 'smooth muscle cell'
        - `BICCN_subclass_label` (8 cell types): 'Oligo', 'Micro-PVM', 'OPC', 'Astro', 'Endo', 'VLMC', 'Peri', 'SMC'
        - `BICCN_cluster_label` (56 cell types)

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered feature_name feature_reference feature_biotype feature_length
feature_id
ENSG00000140463                False         BBS4    NCBITaxon:9606            gene           4148
ENSG00000125533                False      BHLHE23    NCBITaxon:9606            gene           1109
ENSG00000182667                False          NTM    NCBITaxon:9606            gene           6956
ENSG00000101440                False         ASIP    NCBITaxon:9606            gene            845
ENSG00000100568                False        VTI1B    NCBITaxon:9606            gene           5768

adata.var.shape
(14412, 5)
```
- There are 14,412 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata.X.A[1:10, 1:10]
array([[0.       , 1.609438 , 0.       , 0.       , 0.       , 0.       ,
        0.6931472, 0.       , 0.       ],
       [0.       , 1.9459101, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 2.1972246, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 2.0794415, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 2.8903718, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 0.       , 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 2.7080503, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.6931472],
       [0.       , 0.6931472, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ],
       [0.       , 2.4849067, 0.       , 0.       , 0.       , 0.       ,
        0.       , 0.       , 0.       ]], dtype=float32)
# the raw count is in the raw.X layer
adata.raw.X.A[1:10, 1:10]
array([[ 0.,  6.,  0.,  0.,  0.,  0.,  1.,  0.,  0.],
       [ 0.,  3.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0., 10.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  7.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0., 22.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0., 14.,  0.,  0.,  0.,  0.,  0.,  0.,  1.],
       [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],
       [ 0., 11.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]], dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                              nCount_RNA  nFeature_RNA  ...        development_stage observation_joinid
human_AAAGGTATCGGCTGGT-21L8TX_180927_001_A01      5026.0          2228  ...  60-year-old human stage         ikVcksmKO+
human_AAGAACACAACACGAG-21L8TX_180927_001_A01      5798.0          2767  ...  60-year-old human stage         tMqs{UW>Vg
human_AAGCGTTCAAGGAGTC-21L8TX_180927_001_A01      1778.0          1241  ...  60-year-old human stage         @RZf~40Elv
human_AATGGAAAGTCGAGGT-21L8TX_180927_001_A01      4272.0          1990  ...  60-year-old human stage         Wt9Yi`TIFm
human_AATTCCTGTCATCACA-21L8TX_180927_001_A01      3776.0          1720  ...  60-year-old human stage         {tifg;FDSk

adata.obs.shape
(10739, 36)
```
- There are 10,739 cells

**IMPORTANT NOTES 1:** After subsetting for human from the Excitatory 3 species integration and Excitatory 4 species integration, I find that the human cells are the same. Therefore, I will just keep the Excitatroy 3 species integration.

**IMPORTANT NOTES 2:**
- There are 3 files downloaded
```
f8d20729-9133-4113-b49a-3f996e5f8a78.h5ad
62ab3ab7-b0a7-409e-9ed6-4a50d60b964a.h5ad
c0c1469a-1696-450d-949c-7c3a13299ed3.h5ad
```

But the donor id are the same so I'm going to combine the files after computing mean and specifity steps. 
```
import anndata
adata = anndata.read("f8d20729-9133-4113-b49a-3f996e5f8a78.h5ad")
for i in adata.obs['donor_id'].unique():
...     print(i)
...
H18.30.001
H18.30.002
bi005
bi006
F003
F004
F005
F007
F008
M002
M003
M007
M008
F006
M004
M006

#
adata = anndata.read("62ab3ab7-b0a7-409e-9ed6-4a50d60b964a.h5ad")
for i in adata.obs['donor_id'].unique():
...     print(i)
...
H18.30.001
H18.30.002
bi005
bi006
F003
F004
F005
F007
F008
M002
M003
M007
M008
F006
M004
M006

#
adata = anndata.read("c0c1469a-1696-450d-949c-7c3a13299ed3.h5ad")
for i in adata.obs['donor_id'].unique():
...     print(i)
...
H18.30.001
H18.30.002
bi005
bi006
F003
F004
F005
F007
F008
M002
M003
M007
M008
F006
M004
M006
```
