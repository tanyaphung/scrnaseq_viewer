## General information
- Raw counts downloaded on 2023-12-18
```
wget https://datasets.cellxgene.cziscience.com/fafef79b-a863-4e63-a162-dd933d552ace.h5ad
```

## Data exploration
```
adata = anndata.read("fafef79b-a863-4e63-a162-dd933d552ace.h5ad")
```
- The Ensemble gene ID is already included as a row name in `adata.var`
```
>>> adata.var
                 feature_is_filtered feature_name feature_reference feature_biotype feature_length
gene
ENSG00000233576                False      HTR3C2P    NCBITaxon:9606            gene           1057
ENSG00000121410                False         A1BG    NCBITaxon:9606            gene           3999
ENSG00000268895                False     A1BG-AS1    NCBITaxon:9606            gene           3374
ENSG00000148584                False         A1CF    NCBITaxon:9606            gene           9603
ENSG00000175899                False          A2M    NCBITaxon:9606            gene           6318

adata.var.shape
(29813, 5)
```
- There are 29813 genes
- The cell ID is a row name in `adata.obs`
```
adata.obs
                                          Cluster Neighborhood  ...        development_stage observation_joinid
F1S4_160106_003_A01                       Oligo_1         glia  ...        human adult stage         S7hi>`^4kI
F1S4_160106_037_A01                         OPC_2         glia  ...        human adult stage         0USxo>xT5t
F1S4_160106_038_A01                       Astro_5         glia  ...        human adult stage         M;U5_fZ?&N
F1S4_160106_039_A01                       Astro_5         glia  ...        human adult stage         Hk=+!)+&C!
F1S4_160106_040_A01                       Astro_2         glia  ...        human adult stage         7r;h24_bL+

adata.obs.shape
(156285, 27)
```
- There are 156285 cells

1. Need to separate by sequencing technologies

```
adata.obs['assay'].unique()
['Smart-seq v4', '10x 3' v3']
```
2. Cell type annotation
- There seems to be 3 levels, though there does not seem to be a huge difference between the `cell_type` and `Subclass` 
- Level 1: 18 cell types

```
adata.obs['cell_type'].unique()
['oligodendrocyte', 'oligodendrocyte precursor cell', 'astrocyte of the cerebral cortex', 'microglial cell', 'cerebral cortex endothelial cell', ..., 'sncg GABAergic cortical interneuron', 'sst GABAergic cortical interneuron', 'pvalb GABAergic cortical interneuron', 'chandelier pvalb GABAergic cortical interneuron', 'L2/3-6 intratelencephalic projecting glutamat...]
Length: 18
Categories (18, object): ['oligodendrocyte', 'microglial cell', 'oligodendrocyte precursor cell','astrocyte of the cerebral cortex', ..., 'L2/3-6 intratelencephalic projecting glutamat...,'L5 extratelencephalic projecting glutamatergi..., 'vascular leptomeningeal cell','caudal ganglionic eminence derived GABAergic ...]
```
- Level 2: 24 cell types

```
adata.obs['Subclass'].unique()
['Oligo', 'OPC', 'Astro', 'Micro-PVM', 'Endo', ..., 'L5 IT', 'L4 IT', 'L6 IT', 'L6 IT Car3', 'L2/3 IT']
Length: 24
Categories (24, object): ['Astro', 'Chandelier', 'Endo', 'L2/3 IT', ..., 'Sst', 'Sst Chodl', 'VLMC', 'Vip']
```
- Level 3: 151 cell types

```
adata.obs['Cluster'].unique()
['Oligo_1', 'OPC_2', 'Astro_5', 'Astro_2', 'Micro-PVM_2', ..., 'L2/3 IT_2', 'L2/3 IT_9', 'L2/3 IT_7', 'L6 IT Car3_1', 'L5 IT_4']
Length: 151
Categories (151, object): ['Astro_1', 'Astro_2', 'Astro_3', 'Astro_4', ..., 'Vip_20', 'Vip_21',
                           'Vip_22', 'Vip_23']
```