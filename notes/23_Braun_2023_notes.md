## General information
- Paper: https://www.science.org/doi/10.1126/science.adf1226
- Raw counts downloaded on 2025-02-10
```
#ClassOligo
wget https://datasets.cellxgene.cziscience.com/648df3ef-670e-427a-b7b0-305fe05e64c3.h5ad
mv 648df3ef-670e-427a-b7b0-305fe05e64c3.h5ad ClassOligo.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = adata = anndata.read("ClassOligo.h5ad")
```
- Explore `adata.obs`:
```
adata.obs.columns
Index(['Age', 'CellClass', 'Region', 'Subregion', 'organism_ontology_term_id',
       'disease_ontology_term_id', 'self_reported_ethnicity_ontology_term_id',
       'assay_ontology_term_id', 'sex_ontology_term_id',
       'development_stage_ontology_term_id', 'donor_id', 'suspension_type',
       'tissue_type', 'dissection', 'tissue_ontology_term_id',
       'cell_type_ontology_term_id', 'fraction_mitochondrial',
       'fraction_unspliced', 'cell_cycle_score', 'total_genes', 'total_UMIs',
       'sample_id', 'cluster_id', 'is_primary_data', 'cell_type', 'assay',
       'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity',
       'development_stage', 'observation_joinid'],
      dtype='object')
```

- Find content of each column: 
```
for i in adata.obs.columns:
...     print(f"###Column: {i}")
...     print(adata.obs[i].unique())
...
###Column: Age
[8.0, 9.5, 10.0, 8.1, 6.6, ..., 11.5, 13.0, 14.0, 5.5, 7.0]
Length: 17
Categories (17, float64): [5.5, 6.0, 6.6, 6.7, ..., 11.5, 12.0, 13.0, 14.0]
###Column: CellClass
['Oligo']
Categories (1, object): ['Oligo']
###Column: Region
['Midbrain', 'Telencephalon', 'Forebrain', 'Diencephalon', 'Medulla', 'Pons', 'Cerebellum', 'Hindbrain', 'Brain']
Categories (9, object): ['Brain', 'Cerebellum', 'Diencephalon', 'Forebrain', ..., 'Medulla',
                         'Midbrain', 'Pons', 'Telencephalon']
###Column: Subregion
['Midbrain ventral', 'Cortex', 'Subcortex', 'Forebrain', 'Hypothalamus', ..., 'Cerebellum', 'Hindbrain', 'Hippocampus', 'Midbrain dorsal', 'Brain']
Length: 15
Categories (15, object): ['Brain', 'Cerebellum', 'Cortex', 'Forebrain', ..., 'Pons', 'Striatum',
                          'Subcortex', 'Thalamus']
###Column: organism_ontology_term_id
['NCBITaxon:9606']
Categories (1, object): ['NCBITaxon:9606']
###Column: disease_ontology_term_id
['PATO:0000461']
Categories (1, object): ['PATO:0000461']
###Column: self_reported_ethnicity_ontology_term_id
['unknown']
Categories (1, object): ['unknown']
###Column: assay_ontology_term_id
['EFO:0009899', 'EFO:0009922']
Categories (2, object): ['EFO:0009899', 'EFO:0009922']
###Column: sex_ontology_term_id
['unknown']
Categories (1, object): ['unknown']
###Column: development_stage_ontology_term_id
['HsapDv:0000046', 'HsapDv:0000047', 'HsapDv:0000048', 'HsapDv:0000027', 'HsapDv:0000029', ..., 'HsapDv:0000050', 'HsapDv:0000049', 'HsapDv:0000051', 'HsapDv:0000052', 'HsapDv:0000023']
Length: 11
Categories (11, object): ['HsapDv:0000023', 'HsapDv:0000025', 'HsapDv:0000027', 'HsapDv:0000029', ...,
                          'HsapDv:0000049', 'HsapDv:0000050', 'HsapDv:0000051', 'HsapDv:0000052']
###Column: donor_id
['BRC2006', 'XHU:292', 'XHU:297', 'BRC2057', 'BRC2073', ..., 'XDD:359', 'XDD:385', 'XDD:395', 'XDD:400', 'XDD:398']
Length: 23
Categories (23, object): ['BRC2006', 'BRC2057', 'BRC2061', 'BRC2073', ..., 'XHU:292', 'XHU:297',
                          'XHU:305', 'XHU:307']
###Column: suspension_type
['cell']
Categories (1, object): ['cell']
###Column: tissue_type
['tissue']
Categories (1, object): ['tissue']
###Column: dissection
['Ventral midbrain', 'Forebrain cortex', 'Subcortical forebrain', 'Cortex', 'Forebrain', ..., 'Upper cortex', 'Dorsal midbrain', 'Caudate-Putamen', 'Brain', 'Hindbrain']
Length: 25
Categories (25, object): ['Brain', 'Caudate-Putamen', 'Cerebellum', 'Cortex', ...,
                          'Subcortical forebrain', 'Thalamus', 'Upper cortex', 'Ventral midbrain']
###Column: tissue_ontology_term_id
['UBERON:0002474', 'UBERON:0001890', 'UBERON:0002743', 'UBERON:0001851', 'UBERON:0001898', ..., 'UBERON:0016540', 'UBERON:8440044', 'UBERON:0002314', 'UBERON:0005383', 'UBERON:0000955']
Length: 23
Categories (23, object): ['UBERON:0000454', 'UBERON:0000955', 'UBERON:0000988', 'UBERON:0001851', ...,
                          'UBERON:0005383', 'UBERON:0016540', 'UBERON:8440044', 'UBERON:8440051']
###Column: cell_type_ontology_term_id
['CL:0000128']
Categories (1, object): ['CL:0000128']
###Column: fraction_mitochondrial
[0.02424242 0.01268012 0.00337079 ... 0.03860945 0.03008869 0.06177835]
###Column: fraction_unspliced
[0.12354312 0.12795389 0.16179775 ... 0.42829308 0.51874423 0.2594691 ]
###Column: cell_cycle_score
[0.0013986  0.00230548 0.0011236  ... 0.00215351 0.00157789 0.0148376 ]
###Column: total_genes
[1293 1028  644 ... 5859 5798 5258]
###Column: total_UMIs
[ 2145  1735   890 ... 19503 18379 18534]
###Column: sample_id
['10X89_7', '10X89_8', '10X92_1', '10X92_2', '10X92_3', ..., '10X298_3', '10X298_4', '10X302_1', '10X302_3', '10X302_4']
Length: 193
Categories (193, object): ['10X101_1', '10X101_2', '10X101_3', '10X101_4', ..., '10X92_1', '10X92_2',
                           '10X92_3', '10X92_4']
###Column: cluster_id
[2, 1, 5, 4, 6, 7, 3, 0]
Categories (8, uint64): [0, 1, 2, 3, 4, 5, 6, 7]
###Column: is_primary_data
[False]
###Column: cell_type
['oligodendrocyte']
Categories (1, object): ['oligodendrocyte']
###Column: assay
['10x 3' v2', '10x 3' v3']
Categories (2, object): ['10x 3' v2', '10x 3' v3']
###Column: disease
['normal']
Categories (1, object): ['normal']
###Column: organism
['Homo sapiens']
Categories (1, object): ['Homo sapiens']
###Column: sex
['unknown']
Categories (1, object): ['unknown']
###Column: tissue
['cerebellar peduncular complex', 'forebrain', 'basal forebrain', 'cortex', 'hypothalamus', ..., 'occipital cortex', 'upper layers of the cortex', 'midbrain tectum', 'caudate-putamen', 'brain']
Length: 23
Categories (23, object): ['cerebral subcortex', 'brain', 'pons', 'cortex', ..., 'caudate-putamen',
                          'occipital cortex', 'upper layers of the cortex',
                          'lower layers of the cortex']
###Column: self_reported_ethnicity
['unknown']
Categories (1, object): ['unknown']
###Column: development_stage
['9th week post-fertilization stage', '10th week post-fertilization stage', '11th week post-fertilization stage', 'Carnegie stage 20', 'Carnegie stage 22', ..., '13th week post-fertilization stage', '12th week post-fertilization stage', '14th week post-fertilization stage', '15th week post-fertilization stage', 'Carnegie stage 16']
Length: 11
Categories (11, object): ['Carnegie stage 16', 'Carnegie stage 18', 'Carnegie stage 20',
                          'Carnegie stage 22', ..., '12th week post-fertilization stage',
                          '13th week post-fertilization stage', '14th week post-fertilization stage',
                          '15th week post-fertilization stage']
###Column: observation_joinid
['p%u$=bD!E9' 'L=L#e(7hIi' 'd9OIz`_=|r' ... 'i#+jq`sX^g' 'cFwe%Xt#g|'
 '_Z}Jt18-D^']
```

adata.var.head()
                Chromosome    End GeneTotalUMIs  ... feature_biotype feature_length                        feature_type
gene_ids                                         ...
ENSG00000223972       chr1  14409             0  ...            gene            632  transcribed_unprocessed_pseudogene
ENSG00000227232       chr1  29570           263  ...            gene           1351              unprocessed_pseudogene
ENSG00000278267       chr1  17436             0  ...            gene             68                               miRNA
ENSG00000284332       chr1  30503             0  ...            gene            138                               miRNA
ENSG00000268020       chr1  53312             0  ...            gene            840              unprocessed_pseudogene

adata.var['feature_name']
gene_ids
ENSG00000223972                      DDX11L1
ENSG00000227232                       WASH7P
ENSG00000278267                    MIR6859-1
ENSG00000284332                    MIR1302-2
ENSG00000268020                       OR4G4P
                             ...
ENSG00000275757    5_8S_rRNA_ENSG00000275757
ENSG00000278573            ENSG00000278573.1
ENSG00000276017            ENSG00000276017.1
ENSG00000278625           U6_ENSG00000278625
ENSG00000277374           U1_ENSG00000277374
Name: feature_name, Length: 59229, dtype: category
Categories (59229, object): ['5S_rRNA_ENSG00000276861', '5S_rRNA_ENSG00000277411', '5S_rRNA_ENSG00000277488',
                             '5S_rRNA_ENSG00000278457', ..., 'hsa-mir-1253', 'hsa-mir-423', 'hsa-mir-8069-1',
                             'snoZ196']

adata.var.shape
(59229, 11)

```
- There are 59,229 genes

- Count data is store under the X layer
```
adata.X.A
array([[0., 0., 0., ..., 0., 0., 0.],
       [0., 0., 0., ..., 0., 0., 0.],
       [0., 0., 0., ..., 0., 0., 0.],
       ...,
       [0., 1., 0., ..., 0., 0., 0.],
       [0., 0., 0., ..., 0., 0., 0.],
       [0., 0., 0., ..., 0., 0., 0.]], dtype=float32)
```

- The CellID is a row name in `adata.obs`
```
adata.obs.head()
                          Age CellClass  ...                   development_stage observation_joinid
CellID                                   ...
10X89_7:GCGCAGTGTTAAAGAC  8.0     Oligo  ...   9th week post-fertilization stage         p%u$=bD!E9
10X89_7:GTTCTCGCATCACAAC  8.0     Oligo  ...   9th week post-fertilization stage         L=L#e(7hIi
10X89_8:ACGAGCCCAGTACACT  8.0     Oligo  ...   9th week post-fertilization stage         d9OIz`_=|r
10X92_1:AGTGTCACACTTAAGC  9.5     Oligo  ...  10th week post-fertilization stage         2M4~S5@?Hg
10X92_1:CAACCAACACGGTAGA  9.5     Oligo  ...  10th week post-fertilization stage         QmCGKUR9wi

adata.obs.shape
(6190, 33)
```
- There are 6190 cells

## Strategy: 
- Because the All file has >1.5M cells, it's going to be computationally challenging to process this file. Therefore, I am going to download each of the class file
- For each class file (for example ClassNeuron), 
    - separate out for each tissue/region (use keywork `Region`) and age (use keyword `development_stage`)
- For each region and age, concat the h5ad file

## Documents of steps
- Download
```
#ClassNeuron
wget https://datasets.cellxgene.cziscience.com/77bd9788-8b48-4041-b4fd-415b8ad6cf65.h5ad
mv 77bd9788-8b48-4041-b4fd-415b8ad6cf65.h5ad ClassNeuron.h5ad

#ClassRadialGlia
wget https://datasets.cellxgene.cziscience.com/0768590e-3606-40c1-a823-b4fdb89c9f87.h5ad
mv 0768590e-3606-40c1-a823-b4fdb89c9f87.h5ad ClassRadialGlia.h5ad

#ClassNeuroblast
wget https://datasets.cellxgene.cziscience.com/9bdeb1c1-37a9-48c9-8334-a6efebb2eed9.h5ad
mv 9bdeb1c1-37a9-48c9-8334-a6efebb2eed9.h5ad ClassNeuroblast.h5ad

#ClassGlioblastOligo
wget https://datasets.cellxgene.cziscience.com/a517ac77-a03f-4fbe-9c7c-74487fdf75c6.h5ad
mv a517ac77-a03f-4fbe-9c7c-74487fdf75c6.h5ad ClassGlioblastOligo.h5ad

#ClassGlioblast
wget https://datasets.cellxgene.cziscience.com/c300fb0b-3978-481c-ac7c-0e82ac2e2f76.h5ad
mv c300fb0b-3978-481c-ac7c-0e82ac2e2f76.h5ad ClassGlioblast.h5ad

#ClassNeuronalIPC
wget https://datasets.cellxgene.cziscience.com/e8a140cb-6e51-4684-ac67-d60a924c11a6.h5ad
mv e8a140cb-6e51-4684-ac67-d60a924c11a6.h5ad ClassNeuronalIPC.h5ad

#ClassFibroblast
wget https://datasets.cellxgene.cziscience.com/a6827613-dd1f-4811-a823-d920b972fc86.h5ad
mv a6827613-dd1f-4811-a823-d920b972fc86.h5ad ClassFibroblast.h5ad

#ClassErythrocyte
wget https://datasets.cellxgene.cziscience.com/71613e60-9487-4ff1-bd76-cbb00026ef93.h5ad
mv 71613e60-9487-4ff1-bd76-cbb00026ef93.h5ad ClassErythrocyte.h5ad

#ClassImmune
wget https://datasets.cellxgene.cziscience.com/d7a55d23-ca07-446a-9c48-267fe490d1fc.h5ad
mv d7a55d23-ca07-446a-9c48-267fe490d1fc.h5ad ClassImmune.h5ad

#ClassOligo
wget https://datasets.cellxgene.cziscience.com/648df3ef-670e-427a-b7b0-305fe05e64c3.h5ad
mv 648df3ef-670e-427a-b7b0-305fe05e64c3.h5ad ClassOligo.h5ad

#ClassPlacodes
wget https://datasets.cellxgene.cziscience.com/79a6b1ee-309a-4ca6-a8fb-39d2da0e1eb3.h5ad
mv 79a6b1ee-309a-4ca6-a8fb-39d2da0e1eb3.h5ad ClassPlacodes.h5ad

#ClassNeuralCrest
wget https://datasets.cellxgene.cziscience.com/3f9af1de-2ca5-495d-920a-9ac5c903e637.h5ad
mv 3f9af1de-2ca5-495d-920a-9ac5c903e637.h5ad ClassNeuralCrest.h5ad
```

- separate out by tissue and age
```
for type in ClassGlioblast ClassGlioblastOligo ClassRadialGlia ClassNeuron ClassNeuroblast ClassErythrocyte ClassImmune ClassFibroblast ClassOligo ClassNeuronalIPC ClassPlacodes ClassNeuralCrest; do

python split_by_tissue_age.py ${type}.h5ad ${type} summary.txt Preprocessing_scRNA/data/Braun2023_FirstTrimesterBrain/

done;
```

- For each region and age, concat the h5ad file
```
python concat_perTissue_perDevAge.py
```

- Examine `Cerebellum_9wpc.h5ad`
- Make the adata_var file

```
>>> adata = anndata.read("ClassNeuralCrest_Diencephalon_15wpc.h5ad")
>>> adata.var
                       Chromosome     End GeneTotalUMIs  ... feature_biotype feature_length                        feature_type
gene_ids                                                 ...

ENSG00000223972              chr1   14409             0  ...            gene            632  transcribed_unprocessed_pseudogene
ENSG00000227232              chr1   29570            99  ...            gene           1351              unprocessed_pseudogene
ENSG00000278267              chr1   17436             0  ...            gene             68                               miRNA
ENSG00000284332              chr1   30503             0  ...            gene            138                               miRNA
ENSG00000268020              chr1   53312             0  ...            gene            840              unprocessed_pseudogene
...                           ...     ...           ...  ...             ...            ...
...
ENSG00000275757  chr22_KI270733v1  174108             0  ...            gene            153                                rRNA
ENSG00000278573  chr22_KI270734v1   60316             0  ...            gene            603                          pseudogene
ENSG00000276017  chr22_KI270734v1   74814             0  ...            gene           2404                      protein_coding
ENSG00000278625  chrUn_KI270744v1   51114             0  ...            gene            106                               snRNA
ENSG00000277374  chrUn_KI270750v1  148843             0  ...            gene            176                               snRNA

[59229 rows x 11 columns]
>>> adata_var = adata.var[['feature_name']]
>>> adata_var.head()
                feature_name
gene_ids
ENSG00000223972      DDX11L1
ENSG00000227232       WASH7P
ENSG00000278267    MIR6859-1
ENSG00000284332    MIR1302-2
ENSG00000268020       OR4G4P
>>> adata_var.reset_index(inplace=True)
>>> adata_var.head()
          gene_ids feature_name
0  ENSG00000223972      DDX11L1
1  ENSG00000227232       WASH7P
2  ENSG00000278267    MIR6859-1
3  ENSG00000284332    MIR1302-2
4  ENSG00000268020       OR4G4P
>>> adata_var.to_csv("adata_var.txt", sep="\t", index=False)

head adata_var.txt
gene_ids        feature_name
ENSG00000223972 DDX11L1
ENSG00000227232 WASH7P
ENSG00000278267 MIR6859-1
ENSG00000284332 MIR1302-2
ENSG00000268020 OR4G4P
ENSG00000240361 OR4G11P
ENSG00000233750 CICP27
ENSG00000268903 ENSG00000268903.1
ENSG00000269981 ENSG00000269981.1
```

- QC: 
```
python qc_scrna_23_Braun_2023.py 443_Braun2023_Human_2023_FirstTrimester_Brain_CarnegieStage18 Brain CarnegieStage18
```