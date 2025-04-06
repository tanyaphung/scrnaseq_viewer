## General information
- Paper: Transcriptomic cytoarchitecture reveals principles of human neocortex organization
- Link: https://pubmed.ncbi.nlm.nih.gov/37824655/
- What's available to download?
    - For assay 10x 3' v3, there are 8 files corresponding to each of the 8 neocortical regions
    - For assay Smart-seq v4, there is 1 file for 6 tissues
- Raw counts downloaded on 2024-01-10
```
## For 10x
wget https://datasets.cellxgene.cziscience.com/5b347baa-63d2-4b3e-9d17-5d51c0a3dd23.h5ad
mv 5b347baa-63d2-4b3e-9d17-5d51c0a3dd23.h5ad 321_Jorstad2023_V1_Human_2023_10x/original_321_Jorstad2023_V1_Human_2023_10x.h5ad

wget https://datasets.cellxgene.cziscience.com/f4803445-079b-4fbb-b262-ca36c7b6562b.h5ad
mv f4803445-079b-4fbb-b262-ca36c7b6562b.h5ad original_318_Jorstad2023_M1_Human_2023_10x.h5ad

wget https://datasets.cellxgene.cziscience.com/173891cd-830e-4733-861e-5e8ba59973dd.h5ad
mv f4803445-079b-4fbb-b262-ca36c7b6562b.h5ad original_319_Jorstad2023_S1_Human_2023_10x.h5ad

wget https://datasets.cellxgene.cziscience.com/4184721b-0a1b-4bb6-b0fb-dc0bf023acb9.h5ad
mv 4184721b-0a1b-4bb6-b0fb-dc0bf023acb9.h5ad original_320_Jorstad2023_A1_Human_2023_10x.h5ad

wget https://datasets.cellxgene.cziscience.com/26b860c7-77e7-4761-81d5-6ee38bb88117.h5ad
mv 26b860c7-77e7-4761-81d5-6ee38bb88117.h5ad original_322_Jorstad2023_DFC_Human_2023_10x.h5ad

wget https://datasets.cellxgene.cziscience.com/5b37a4b8-7bbd-456c-a0ba-40c78411a608.h5ad
mv 5b37a4b8-7bbd-456c-a0ba-40c78411a608.h5ad original_323_Jorstad2023_ACC_Human_2023_10x.h5ad

wget https://datasets.cellxgene.cziscience.com/fffa9200-0e57-4f8c-b649-d4602434650b.h5ad
mv fffa9200-0e57-4f8c-b649-d4602434650b.h5ad original_324_Jorstad2023_MTG_Human_2023_10x.h5ad

wget https://datasets.cellxgene.cziscience.com/36f4502e-4566-45aa-a20c-cfdc0714f601.h5ad
mv 36f4502e-4566-45aa-a20c-cfdc0714f601.h5ad original_325_Jorstad2023_AnG_Human_2023_10x.h5ad

## For SMART-seq
wget https://datasets.cellxgene.cziscience.com/bae0de2c-d2d9-4dd3-ac78-49ca57f20929.h5ad
mv bae0de2c-d2d9-4dd3-ac78-49ca57f20929.h5ad Jorstadetal2023_Neocortex_smartseq.h5ad
```

## Data exploration
### 10x
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("5b347baa-63d2-4b3e-9d17-5d51c0a3dd23.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'Class', 'CrossArea_subclass', 'CrossArea_cluster',
       'WithinArea_subclass', 'WithinArea_cluster', 'Location', 'Region', 'Subregion',  'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'development_stage'
    - Location: caudal medial
    - Region: V1 (primary visual cortex)
    - Subregion: V1C, V1_L5, V1C_l5
    - Normal (no disease)
    - Assay (10x 3' v3)
    - Sex: male and female
    - Organism: Homo sapiens
    - Tissue: primary visual cortex
    - Development stage: adults: '50-year-old human stage', '29-year-old human stage', '42-year-old human stage', '43-year-old human stage', '60-year-old human stage'
    - How many cell types? 
        - `Class` (3 cell types): 'excitatory', 'inhibitory', 'non-neuronal' -> level_1
        - `CrossArea_subclass` (24 cell types): 'Lamp5', 'Sncg', 'Lamp5 Lhx6', 'Pax6', 'Vip', 'L5 ET', 'L5/6 NP', 'L6 CT', 'L6b', 'Astro', 'Endo', 'OPC', 'VLMC', 'Micro/PVM', 'Oligo', 'L2/3 IT', 'L4 IT', 'L6 IT Car3', 'L6 IT', 'L5 IT', 'Chandelier', 'Pvalb', 'Sst', 'Sst Chodl'
        - `CrossArea_cluster` (146 cell types)
        - `WithinArea_subclass` (24 cell ytpes): same as `CrossArea_subclass` -> level_2
        - `WithinArea_cluster` (130 cell types) -> level_3
        - `cell_type` (18 cell types): 'lamp5 GABAergic cortical interneuron', 'sncg GABAergic cortical interneuron', 'caudal ganglionic eminence derived GABAergic cortical interneuron', 'vip GABAergic cortical interneuron', 'L5 extratelencephalic projecting glutamatergic cortical neuron', 'near-projecting glutamatergic cortical neuron', 'corticothalamic-projecting glutamatergic cortical neuron', 'L6b glutamatergic cortical neuron', 'astrocyte of the cerebral cortex', 'cerebral cortex endothelial cell', 'oligodendrocyte precursor cell', 'vascular leptomeningeal cell', 'microglial cell', 'oligodendrocyte', 'L2/3-6 intratelencephalic projecting glutamatergic cortical neuron', 'chandelier pvalb GABAergic cortical interneuron', 'pvalb GABAergic cortical interneuron', 'sst GABAergic cortical interneuron'
        - 

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered feature_name feature_reference
gene
ENSG00000233576                False      HTR3C2P    NCBITaxon:9606  \
ENSG00000121410                False         A1BG    NCBITaxon:9606
ENSG00000268895                False     A1BG-AS1    NCBITaxon:9606
ENSG00000148584                False         A1CF    NCBITaxon:9606
ENSG00000175899                False          A2M    NCBITaxon:9606

                feature_biotype feature_length
gene
ENSG00000233576            gene           1057
ENSG00000121410            gene           3999
ENSG00000268895            gene           3374
ENSG00000148584            gene           9603
ENSG00000175899            gene           6318

adata.var.shape
(29444, 5)
```
- There are 29,444 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata.X.A[1:10, 1:10]
[[0.        0.        0.        0.        0.        0.        0.
  0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  3.6604488 0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.       ]
 [0.        4.3441725 0.        0.        0.        0.        0.
  0.        0.       ]
 [4.1204624 0.        0.        0.        0.        4.1204624 0.
  0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.       ]
 [4.298724  0.        0.        0.        0.        0.        0.
  0.        0.       ]]
# the raw count is in the raw.X layer
adata.raw.X.A[1:10, 1:10]
[[0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 1. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 1. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0.]]
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                           Class CrossArea_subclass  ...        development_stage observation_joinid
AACGGGACACCATTCC-9L8TX_200107_01_C10  inhibitory              Lamp5  ...  50-year-old human stage         k`A9@Jt{9=
ACTCTCGAGCTACAAA-9L8TX_200107_01_C10  inhibitory              Lamp5  ...  50-year-old human stage         3lORk`ETj(
AGAGCAGAGTGGTTGG-9L8TX_200107_01_C10  inhibitory              Lamp5  ...  50-year-old human stage         J_|_BAeo!B
AGGTCTAGTGGACCTC-9L8TX_200107_01_C10  inhibitory              Lamp5  ...  50-year-old human stage         qziY9&X)Og
AGGTTGTTCCGAGTGC-9L8TX_200107_01_C10  inhibitory              Lamp5  ...  50-year-old human stage         MDg5VsRjEQ

[5 rows x 33 columns]

adata.obs.shape
(241077, 33)
```
- There are 241,077 cells

### SMART-seq
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("bae0de2c-d2d9-4dd3-ac78-49ca57f20929.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'cluster', 'class', 'subclass', 'region', 'cortical_layer', 'cause_of_death', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'development_stage', 'observation_joinid'
    - region: (8) 'A1C', 'CgG', 'M1lm', 'M1ul', 'MTG', 'S1lm', 'S1ul', 'V1C'
    - cortical_layer: (13) 
    - cause of death: 'cadiovascular disease', 'mitral valve prolapse'
    - assay: Smart-seq v4
    - disease: normal
    - sex: female and male
    - tissue: (6) 'primary motor cortex', 'primary visual cortex', 'anterior cingulate gyrus',
                         'middle temporal gyrus', 'primary somatosensory cortex',
                         'primary auditory cortex'
    - development stage: adults: '43-year-old human stage', '50-year-old human stage', '54-year-old human stage'
    - How many cell types? 
        - `class` (3 cell types)
        - `subclass` (19 cell types)
        - `cluster` (121 cell types)
        - `cell_type` (19 cell types)

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered feature_name feature_reference feature_biotype feature_length
ENSG00000233576                False      HTR3C2P    NCBITaxon:9606            gene           1057
ENSG00000121410                False         A1BG    NCBITaxon:9606            gene           3999
ENSG00000268895                False     A1BG-AS1    NCBITaxon:9606            gene           3374
ENSG00000148584                False         A1CF    NCBITaxon:9606            gene           9603
ENSG00000175899                False          A2M    NCBITaxon:9606            gene           6318

adata.var.shape
(35812, 5)
```
- There are 35,812 genes

- Count data is store under the X layer:

```
# the X layer shows integers
adata.X.A[1:10, 1:10]
array([[118.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.],
       [  0.,   0.,   0.,   0.,   0.,   6.,   0.,   0.,   0.],
       [159.,   0.,   0.,   0.,   0.,   2.,   0.,   0.,   1.],
       [  1.,   0.,   0.,   1.,   0.,   1.,   0.,   0.,   1.],
       [ 10.,   0.,   0.,   0.,   0.,   1.,   0.,   3.,   1.],
       [  0.,   0.,   0., 231.,   0.,   1.,   0.,   0.,   0.],
       [  1.,   0.,   0.,   0.,   0.,   1.,   0.,   0.,   0.],
       [  0.,   0.,   0.,   0.,   0.,   0.,   0.,   1.,   1.],
       [  0.,   0.,   0.,   0.,   0.,   1.,   0.,   0.,   0.]],
      dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                    suspension_type                cluster      class  ... self_reported_ethnicity        development_stage observation_joinid
F2S4_160113_027_A01         nucleus                    NaN        NaN  ...                European  50-year-old human stage         $j-^M1E9?*
F2S4_160113_027_B01         nucleus      Inh L2-5 VIP TOX2  GABAergic  ...                European  50-year-old human stage         g~daSNKdwS
F2S4_160113_027_C01         nucleus     Inh L1 LAMP5 GGT8P  GABAergic  ...                European  50-year-old human stage         ;{0&$@~aP~
F2S4_160113_027_D01         nucleus      Inh L1 LAMP5 NDNF  GABAergic  ...                European  50-year-old human stage         ?KjSxtN%(K
F2S4_160113_027_E01         nucleus  Inh L1-3 VIP ZNF322P1  GABAergic  ...                European  50-year-old human stage         LE{T(vw2u~

adata.obs.shape
(49417, 38)
```
- There are 49,417 cells