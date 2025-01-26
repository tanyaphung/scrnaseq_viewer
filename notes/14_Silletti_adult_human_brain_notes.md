## General information
- Example for 7_Siletti_CerebralCortex.PrCG.M1C_Human_2022
- Raw counts downloaded on 2023-04-05
```
curl -o local.h5ad "https://corpora-data-prod.s3.amazonaws.com/e3ff5787-ab6b-49fe-9e26-68e4b3904042/local.h5ad?AWSAccessKeyId=ASIATLYQ5N5XRDL55NOZ&Signature=Jn5Qm9O570w3%2FMqC3GNcsrLrCUE%3D&x-amz-security-token=IQoJb3JpZ2luX2VjEOz%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLXdlc3QtMiJHMEUCIQC%2FeF5Wq%2BBJmW03qHEtHA9245envwex7vaaGDcj3IjXBwIgSgYobQuQDmbjRCbctDUUmlkSYEoseURZ7fj9OS092Asq9AMI9f%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FARABGgwyMzE0MjY4NDY1NzUiDLWHUHK%2BX5xivjDkNSrIAxfNmQIEvMuNjTz4Qjfq1bgAFSB6GQ6g6%2F5%2FMWs6%2BlCI%2FYs2FwzqLi6tpXSd7IJLPugNRTOqoS3Oq6IaNDQbmSrHLxAUTw0hIYPvV9pp2%2BwSOhgZyovCNTaxaW8cNMdwFTbXgtjrbx4eL%2BdDyolQZUd6vd2B8Lqyqtdqn5fqzW9kfWE87hHtLaB9J0L%2Fp%2By9V2HQ1x4zncrToVlPcHTNY%2B6aRZKHmOLAg439MjO4R7MzFuzEzDk%2BZYXP%2BRYf9k19U%2FrpNxqqtgzmqq6JsL88Sb1YL7K1kX1yQrUzEW0bPlzciMB7VUxYjvoqMS0lIYRZhtt0eVfix0KREzAyHboq3n2pOjHwOQ1qv7m2U%2F6pwF1vaeDSeO1aBTB8k8nleeF9WJ16JcZyuAlJzVjVqjfrjjkhWRyZcThjsOE%2B6K6QVXUmsceW%2F05ohj%2BhqizC2Z%2BuDzSuoEA4L1jHFXwis%2BBNitg9Z4mUCVl5%2B%2BMd5jZPizRjdp4N5bUodTvbgWRVxF5V6vAiVF9zZHOGQULZMWq5OdIYDV4%2FXamGstpYLlXLENli%2FNS2izU%2BRetrnHpmk62qCxYsOE3vgEwqfMBqjkcxcZLkoVfes9CebTDTmNCiBjqlAcysViDhajf2ySanZxi3h3kPiNrxubKnJWee8heMKZXyO2lHHZT4IMdmRKtcQfiwjeTHRG%2Bw%2BjDBNA8pur6F0dsk%2BivamXrTsSMrcndHdpnG%2BbVY7jvHlzhzoI%2BrHo7KqLfKeRmXQh7FJ%2BmKgYuQCT8KBx0czBCv3AiIi0SZL%2Farnjt3X4w90wf3lsDw89TTJH4EYKywVkvHGjoN4tRboUC7fUcvPA%3D%3D&Expires=1683838735"
```
- The downloaded data is called `local.h5ad`

## Data exploration

```
adata = anndata.read("local.h5ad")
```

- Explore `adata.obs`: 'dissection', 'supercluster_term', 'cell_type', 'assay', 'disease', 'organism', 'sex',
       'tissue', 'self_reported_ethnicity', 'development_stage'
    - dissection: 'Hypothalamus (HTH) - supraoptic region of HTH...
    - assay: 10x 3' v3
    - disease: normal
    - organism: 'Homo sapiens'
    - sex: male
    - tissue: hypothalamus
    - development_stage: '50-year-old human stage'
    - How many cell types? 
        - `cell_type` (11 cell types): 'neuron', 'oligodendrocyte', 'oligodendrocyte precursor cell', 'ependymal cell', 'astrocyte', ..., 'pericyte', 'vascular associated smooth muscle cell', 'endothelial cell', 'central nervous system macrophage', 'leukocyte'
        - `supercluster_term` (21 cell types): 'Upper-layer intratelencephalic', 'Hippocampal CA1-3', 'Amygdala excitatory', 'Miscellaneous', 'Deep-layer intratelencephalic', 'Deep-layer near-projecting', 'Deep-layer corticothalamic and 6b', 'Splatter', 'Eccentric medium spiny neuron', 'Medium spiny neuron', 'MGE interneuron', 'LAMP5-LHX6 and Chandelier', 'CGE interneuron', 'Committed oligodendrocyte precursor', 'Oligodendrocyte', 'Oligodendrocyte precursor', 'Ependymal', 'Astrocyte', 'Fibroblast', 'Vascular', 'Microglia'
- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                Biotype Chromosome        End  ...   feature_name feature_reference  feature_biotype
ensembl_ids                                    ...
ENSG00000135406     n/a      chr12   49298686  ...           PRPH    NCBITaxon:9606             gene
ENSG00000139211     n/a      chr12   47079951  ...         AMIGO2    NCBITaxon:9606             gene
ENSG00000164616     n/a       chr5  135951591  ...        FBXL21P    NCBITaxon:9606             gene
ENSG00000234715     n/a       chr7  103514007  ...   CTB-107G13.1    NCBITaxon:9606             gene
ENSG00000232896     n/a       chr9  110605219  ...  RP11-410K21.2    NCBITaxon:9606             gene

[5 rows x 9 columns]

adata.var.shape
(59357, 9)
```
- There are 59,357 genes

- Count data is store under the X layer:

```
# the X layer shows raw count
adata.X.A[1:15, 1:15]
array([[0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 5., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 8., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 1., 0., 0., 2., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2., 0.],
       [0., 0., 1., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
       [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]],
      dtype=float32)
```

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                 roi organism_ontology_term_id  ... self_reported_ethnicity        development_stage
CellID                                                          ...

10X270_3:CCACTTGGTCTACAAC  Human M1C            NCBITaxon:9606  ...                European  42-year-old human stage
10X159_2:GCTCAAAAGTAGTCAA  Human M1C            NCBITaxon:9606  ...                European  60-year-old human stage
10X159_1:CCTCTCCTCCGTATAG  Human M1C            NCBITaxon:9606  ...                European  60-year-old human stage
10X159_3:GGTTAACTCTCAACCC  Human M1C            NCBITaxon:9606  ...                European  60-year-old human stage
10X160_8:AGCGATTAGCAAGTCG  Human M1C            NCBITaxon:9606  ...                European  60-year-old human stage

[5 rows x 30 columns]

adata.obs.shape
(116576, 30)
```
- There are 116576 cells

## Preprocessing
- `code/preprocessing/14_Siletti_AdultHumanBrain/qc_scrna_siletti_adulthumanbrain.py`
- I used a snakemake pipeline to run this python script in parallel: `code/preprocessing/14_Siletti_AdultHumanBrain/snakefile_qc_scrna_siletti.smk`
- update cell type columns (this is because when this data was originally processed, we did not come to the conclusion in our group that we would rename the first level of cell type annotation to cell_type_level_1)
- `code/preprocessing/14_Siletti_AdultHumanBrain/update_celltype.py` (run with `code/preprocessing/14_Siletti_AdultHumanBrain/run_update_celltype.sh`)

## Processing for FUMA celltype
### Mean gene expression per cell type
- Python script: `code/postprocessing/compute_sumstat_magma.py`
       - I  used a snakemake pipeline to run the above Pythons script: `code/postprocessing/compute_sumstat_magma.smk`