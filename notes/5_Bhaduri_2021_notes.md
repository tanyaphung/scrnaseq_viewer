## General information
- Paper: An atlas of cortical arealization identifies dynamic molecular signatures
- Link: https://www.nature.com/articles/s41586-021-03910-8#Abs1
- Raw counts downloaded on 2024-01-08
```
wget https://datasets.cellxgene.cziscience.com/d0bbb292-7b5d-4f49-b3d8-046bc32a444c.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("d0bbb292-7b5d-4f49-b3d8-046bc32a444c.h5ad")
```
- Explore `adata.obs`:
    - Relevant columns: 'donor_id', 'brain_region', 'cortical_area', 'cluster_label', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'development_stage'
    - Normal (no disease)
    - Assay (10x 3' v2)
    - Sex: unknown
    - Organism: Homo sapiens
    - brain_region: neocortex
        - 5 different cortical areas: 'pfc', 'motor', 'parietal', 'somatosensory', 'temporal', 'v1'. However this is also annotated in the column `tissue` (but there are 6): 'prefrontal cortex', 'primary motor cortex', 'parietal cortex', 'primary somatosensory cortex', 'temporal cortex', 'primary visual cortex'
    - Development stage: '14th week post-fertilization human stage', '16th week post-fertilization human stage', '17th week post-fertilization human stage', '18th week post-fertilization human stage', '19th week post-fertilization human stage', '20th week post-fertilization human stage', '22nd week post-fertilization human stage', '25th week post-fertilization human stage'
    - How many cell types? There are 2 colnames for cell type: 
        - `cell_type` (10 cell type): 'progenitor cell', 'cerebral cortex GABAergic interneuron', 'native cell', 'oligodendrocyte precursor cell', 'forebrain radial glial cell', 'microglial cell', 'blood vessel endothelial cell', 'Cajal-Retzius cell', 'glutamatergic neuron', 'cerebral cortex endothelial cell'
        - `cluster_label` (209 cell types)

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                 feature_is_filtered  feature_name feature_reference
ENSG00000238009                False  RP11-34P13.7    NCBITaxon:9606  \
ENSG00000228463                False    AP006222.1    NCBITaxon:9606
ENSG00000237094                False  RP4-669L17.4    NCBITaxon:9606
ENSG00000237491                False     LINC01409    NCBITaxon:9606
ENSG00000225880                False     LINC00115    NCBITaxon:9606

                feature_biotype feature_length
ENSG00000238009            gene           3726
ENSG00000228463            gene           8224
ENSG00000237094            gene           6204
ENSG00000237491            gene           8413
ENSG00000225880            gene           1317

adata.var.shape
(26902, 5)
```
- There are 26,902 genes

- Count data is store under the raw.X layer:

```
# the X layer shows decimals
adata.X[1:10, 1:10]
[[-0.319264   -0.07110304 -0.12303548 -0.0975266  -0.08828668 -0.02314478
  -0.08136759 -0.01465801 -0.03117384]
 [-0.319264   -0.07110304 -0.12303548 -0.0975266  -0.08828668 -0.02314478
  -0.08136759 -0.01465801 -0.03117384]
 [ 1.2796279  -0.07110304 -0.12303548 -0.0975266  -0.08828668 -0.02314478
  -0.08136759 -0.01465801 -0.03117384]
 [-0.319264   -0.07110304 -0.12303548 -0.0975266  -0.08828668 -0.02314478
  -0.08136759 -0.01465801 -0.03117384]
 [-0.319264   -0.07110304 -0.12303548 -0.0975266  -0.08828668 -0.02314478
  -0.08136759 -0.01465801 -0.03117384]
 [ 1.0217204  -0.07110304 -0.12303548 -0.0975266  -0.08828668 -0.02314478
  -0.08136759 -0.01465801 -0.03117384]
 [-0.319264   -0.07110304 -0.12303548 -0.0975266  -0.08828668 -0.02314478
  -0.08136759 -0.01465801 -0.03117384]
 [-0.319264   -0.07110304 -0.12303548 -0.0975266  -0.08828668 -0.02314478
  -0.08136759 -0.01465801 -0.03117384]
 [-0.319264   -0.07110304 -0.12303548 -0.0975266  -0.08828668 -0.02314478
  -0.08136759 -0.01465801 -0.03117384]]
# the raw count is in the raw.X layer
adata.raw.X.A[1:10, 1:10]
[[0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]]
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                            development_stage_ontology_term_id donor_id
gw19_2_PFC_AAACCCAGTACAAACA                     HsapDv:0000056   GW19_2  \
gw19_2_PFC_AAACCCAGTCTTGCGG                     HsapDv:0000056   GW19_2
gw19_2_PFC_AAACGAAAGCACCAGA                     HsapDv:0000056   GW19_2
gw19_2_PFC_AAACGAAAGTGTACCT                     HsapDv:0000056   GW19_2
gw19_2_PFC_AAACGAACACCCTAAA                     HsapDv:0000056   GW19_2

adata.obs.shape
(457965, 31)
```
- There are 457,965 cells
- NOTE: stratify by cortical area and developmental age