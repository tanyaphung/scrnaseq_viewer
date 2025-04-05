## General information
- In June 2023 I downloaded the reference data from SEA-AD MTG from Allen Brain website. 
- Note that on the cellxgene page (https://cellxgene.cziscience.com/collections/1ca90a2d-2943-483d-b678-b809bf464c30) there are data from the `dorsolateral prefrontal cortex` and `middle temporal gyrus` that includes both dementia and normal. Typically I would process these by subsetting for the normal only. However, in this case, the file is so large (each file has more than 1million cells). Therefore, at this time, I can't process this. Rather, I am using the scRNAseq data that I downloaded from the reference from the Allen Brain Atlas website (http://portal.brain-map.org/atlases-and-data/rnaseq/human-mtg-10x_sea-ad)
- Paper: Integrated multimodal cell atlas of Alzheimer’s disease
- Link: https://www.biorxiv.org/content/10.1101/2023.05.08.539485v2
- Raw counts downloaded on 2023-06-12

### Examine the data

```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("Reference_MTG_RNAseq_all-nuclei.2022-06-07.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'species_label', 'age_label', 'region_label', 'cortical_layer_label', 'cluster_label', 'subclass_label', 'class_label'
        - species_label: 'Homo Sapiens'
        - age_label: ['50 yrs', '42 yrs', '29 yrs', '43 yrs', '60 yrs']
        - region_label: 'MTG'
        - cortical_layer_label: 'all', 'L5'
        - 
    - How many cell types?
        - `class_label` (3 cell types)
        - `subclass_label` (25 cell types)
        - `cluster_label` (128 cell types)

```
adata = anndata.read("Reference_MTG_RNAseq_all-nuclei.2022-06-07.h5ad", backed="r")
<!-- adata
AnnData object with n_obs × n_vars = 166868 × 36601 backed at 'Reference_MTG_RNAseq_all-nuclei.2022-06-07.h5ad'
    obs: 'sample_name', 'donor_sex_label', 'external_donor_name_label', 'species_label', 'age_label', 'region_label', 'cortical_layer_label', 'full_genotype_label', 'QCpass', 'cluster_label', 'cluster_confidence', 'subclass_label', 'subclass_confidence', 'class_label', 'class_confidence', 'GA_QCpass', 'GA_cluster_label', 'GA_subclass_label', 'GA_neighborhood_label', 'CA_QCpass', 'CA_cluster_label', 'CA_subclass_label', 'CA_neighborhood_label', 'cluster_color', 'cluster_order', 'subclass_color', 'subclass_order', 'class_color', 'class_order', 'GA_cluster_color', 'GA_cluster_order', 'GA_subclass_color', 'GA_subclass_order', 'CA_cluster_color', 'CA_cluster_order', 'CA_subclass_color', 'CA_subclass_order', 'cell_type_accession_label' -->
```
- Examine the matrix
```
adata.X[100,1500:2000].A
array([[0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0., 0., 0., 0., 1., 0.,
        0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 2., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        1., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0.,
        1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 1., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 1.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 2., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 2., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 4., 0., 0., 0., 2., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 1., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0., 1., 0., 1.,
        0., 1., 0., 0., 0., 0., 0., 1., 0., 1., 0., 0., 1., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 1., 0., 0., 0., 0.,
        0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 1., 0.,
        0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        0., 0., 0., 0.]], dtype=float32)
```
    + Since the values are integers, I assume that this is the raw count

- Examine the observation
```
adata.obs
                                                                                   sample_name donor_sex_label  ... CA_subclass_order cell_type_accession_label
specimen_name                                                                                                   ...                                          
AAACCCACAACTCATG-LKTX_191204_01_A01-1156636525  AAACCCACAACTCATG-LKTX_191204_01_A01-1156636525               M  ...              18.0             CS202204130_8
AAACCCACAATAGTAG-LKTX_191204_01_A01-1156636525  AAACCCACAATAGTAG-LKTX_191204_01_A01-1156636525               M  ...               NaN                       NaN
AAACCCACACGGTGTC-LKTX_191204_01_A01-1156636525  AAACCCACACGGTGTC-LKTX_191204_01_A01-1156636525               M  ...               7.0           CS202204130_105
AAACCCACACTCTGCT-LKTX_191204_01_A01-1156636525  AAACCCACACTCTGCT-LKTX_191204_01_A01-1156636525               M  ...              15.0            CS202204130_89
AAACCCACATCAGCAT-LKTX_191204_01_A01-1156636525  AAACCCACATCAGCAT-LKTX_191204_01_A01-1156636525               M  ...               9.0            CS202204130_93
...                                                                                        ...             ...  ...               ...                       ...
TTGGGCGTCCTCTTTC-L8TX_200107_01_A09-1156636575  TTGGGCGTCCTCTTTC-L8TX_200107_01_A09-1156636575               F  ...              15.0                       NaN
TTGGGCGTCGTGGGTC-L8TX_200107_01_A09-1156636575  TTGGGCGTCGTGGGTC-L8TX_200107_01_A09-1156636575               F  ...              15.0                       NaN
TTGGGTAGTTTGGGAG-L8TX_200107_01_A09-1156636575  TTGGGTAGTTTGGGAG-L8TX_200107_01_A09-1156636575               F  ...              14.0                       NaN
TTGTGGAAGTGGACGT-L8TX_200107_01_A09-1156636575  TTGTGGAAGTGGACGT-L8TX_200107_01_A09-1156636575               F  ...              14.0                       NaN
TTTCGATAGGAATTAC-L8TX_200107_01_A09-1156636575  TTTCGATAGGAATTAC-L8TX_200107_01_A09-1156636575               F  ...              23.0                       NaN

[166868 rows x 38 columns]

# find columns
adata.obs.columns
Index(['sample_name', 'donor_sex_label', 'external_donor_name_label',
       'species_label', 'age_label', 'region_label', 'cortical_layer_label',
       'full_genotype_label', 'QCpass', 'cluster_label', 'cluster_confidence',
       'subclass_label', 'subclass_confidence', 'class_label',
       'class_confidence', 'GA_QCpass', 'GA_cluster_label',
       'GA_subclass_label', 'GA_neighborhood_label', 'CA_QCpass',
       'CA_cluster_label', 'CA_subclass_label', 'CA_neighborhood_label',
       'cluster_color', 'cluster_order', 'subclass_color', 'subclass_order',
       'class_color', 'class_order', 'GA_cluster_color', 'GA_cluster_order',
       'GA_subclass_color', 'GA_subclass_order', 'CA_cluster_color',
       'CA_cluster_order', 'CA_subclass_color', 'CA_subclass_order',
       'cell_type_accession_label'],
      dtype='object')
```
    + `class_label`: 3 cell types, exluding nan (['Neuronal: GABAergic', nan, 'Neuronal: Glutamatergic', 'Non-neuronal and Non-neural'])
    + `subclass_label`: 24 cell types (excluding nan)
    + `cluster_label`: 127 cell types (excluding nan)

- Examine the var
```
adata.var
Empty DataFrame
Columns: []
Index: [MIR1302-2HG, FAM138A, OR4F5, AL627309.1, AL627309.3, AL627309.2, AL627309.5, AL627309.4, AP006222.2, AL732372.1, OR4F29, AC114498.1, OR4F16, AL669831.2, LINC01409, FAM87B, LINC01128, LINC00115, FAM41C, AL645608.6, AL645608.2, AL645608.4, LINC02593, SAMD11, NOC2L, KLHL17, PLEKHN1, PERM1, AL645608.7, HES4, ISG15, AL645608.1, AGRN, AL645608.5, AL645608.8, RNF223, C1orf159, AL390719.3, LINC01342, AL390719.2, TTLL10-AS1, TTLL10, TNFRSF18, TNFRSF4, SDF4, B3GALT6, C1QTNF12, AL162741.1, UBE2J2, LINC01786, SCNN1D, ACAP3, PUSL1, INTS11, AL139287.1, CPTP, TAS1R3, DVL1, MXRA8, AURKAIP1, CCNL2, MRPL20-AS1, MRPL20, AL391244.2, ANKRD65, AL391244.1, TMEM88B, LINC01770, VWA1, ATAD3C, ATAD3B, ATAD3A, TMEM240, SSU72, AL645728.1, FNDC10, AL691432.4, AL691432.2, MIB2, MMP23B, CDK11B, FO704657.1, SLC35E2B, CDK11A, SLC35E2A, NADK, GNB1, AL109917.1, CALML6, TMEM52, CFAP74, AL391845.2, GABRD, AL391845.1, PRKCZ, AL590822.2, PRKCZ-AS1, FAAP20, AL590822.1, SKI, ...]

[36601 rows x 0 columns]
```
    + This is an empty dataframe with the row index being the gene symbol