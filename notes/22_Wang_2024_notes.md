# General information
- Paper: https://www.biorxiv.org/content/10.1101/2024.01.16.575956v2
- Raw counts downloaded on 2025-02-09

```
wget https://datasets.cellxgene.cziscience.com/23f48c19-6c88-431c-b518-8d8be7ac8706.h5ad
```

## Data exploration
```
import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
adata = anndata.read("23f48c19-6c88-431c-b518-8d8be7ac8706.h5ad")
```
- Explore `adata.obs`:
    - All columns: 'Sample_ID', 'Estimated_postconceptional_age_in_days', 'Group',
       'Region', 'nCount_RNA', 'nFeature_RNA', 'ATAC_fragments_in_peaks',
       'Percentage_reads_in_peaks', 'TSS.enrichment', 'Nucleosome_signal',
       'Scrublet_doublet_score', 'S.Score', 'G2M.Score', 'Class', 'Subclass',
       'Type', 'Type_updated', 'Cluster', 'organism_ontology_term_id',
       'tissue_ontology_term_id', 'tissue_type', 'assay_ontology_term_id',
       'disease_ontology_term_id', 'cell_type_ontology_term_id',
       'self_reported_ethnicity_ontology_term_id',
       'development_stage_ontology_term_id', 'sex_ontology_term_id',
       'donor_id', 'suspension_type', 'is_primary_data', 'cell_type', 'assay',
       'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity',
       'development_stage', 'observation_joinid'

    - Find content of each column: 
    ```
        for i in adata.obs.columns:
        ...     print(i)
        ...     print(adata.obs[i].unique())
        ...
        Sample_ID
        ['ARKFrozen-18-PFC', 'ARKFrozen-19-V1', 'ARKFrozen-1-PFC', 'ARKFrozen-20-CTX', 'ARKFrozen-32-V1', ..., 'NIH-5554-BA17', 'NIH-5554-BA9', 'NIH-5900-BA17', 'NIH-M1154-BA10-2', 'NIH-M2837-BA10-2']
        Length: 38
        Categories (38, object): ['ARKFrozen-1-PFC', 'ARKFrozen-8-V1', 'ARKFrozen-18-PFC', 'ARKFrozen-19-V1',
                                ..., 'NIH-5554-BA17', 'NIH-5900-BA17', 'NIH-M1154-BA10-2',
                                'NIH-M2837-BA10-2']
        Estimated_postconceptional_age_in_days
        [  98  126   84   91  112  147   77  178   60   63   54  455  574  226
        5353  373  507  498 5058 5119 5020  224  239  246]
        Group
        ['Second_trimester', 'Third_trimester', 'First_trimester', 'Infancy', 'Adolescence']
        Categories (5, object): ['Adolescence', 'First_trimester', 'Infancy', 'Second_trimester',
                                'Third_trimester']
        Region
        ['PFC', 'V1', 'General']
        Categories (3, object): ['General', 'PFC', 'V1']
        nCount_RNA
        [ 5071  1106  5564 ... 20523 18427 23880]
        nFeature_RNA
        [2436  786 2508 ... 6766 7297 7014]
        ATAC_fragments_in_peaks
        [ 6329   381  6727 ... 12971 17983 18991]
        Percentage_reads_in_peaks
        [52.29714097 50.19762846 41.99388226 ... 28.45508886 23.18141628
        27.50949275]
        TSS.enrichment
        [ 4.27985299 11.07463965  4.3300329  ...  4.13977167  4.2580707
        3.94519274]
        Nucleosome_signal
        [0.67565217 0.54081633 0.72283951 ... 1.09235936 1.1864813  1.17435321]
        Scrublet_doublet_score
        [0.07869703 0.01824432 0.11840563 ... 0.00350373 0.00320798 0.002917  ]
        S.Score
        [ 0.02841498  0.00689511 -0.08571762 ... -0.00475765 -0.01273367
        0.05508737]
        G2M.Score
        [-0.04040746 -0.04844797 -0.00777326 ... -0.00237759  0.03110667
        -0.06050836]
        Class
        ['Neuron', 'Progenitor', 'Glia', 'Unknown', 'Immune', 'Vascular']
        Categories (6, object): ['Glia', 'Immune', 'Neuron', 'Progenitor', 'Unknown', 'Vascular']
        Subclass
        ['Glutamatergic neuron', 'Radial glia', 'GABAergic neuron', 'IPC-EN', 'Astrocyte', ..., 'Microglia', 'Vascular', 'Oligodendrocyte', 'Cajal-Retzius cell', 'IPC-Glia']
        Length: 12
        Categories (12, object): ['Astrocyte', 'Cajal-Retzius cell', 'GABAergic neuron',
                                'Glutamatergic neuron', ..., 'Oligodendrocyte', 'Radial glia', 'Unknown', 'Vascular']
        Type
        ['EN-IT-Immature', 'EN-Newborn', 'EN-L5-ET', 'RG-vRG', 'IN-MGE-Immature', ..., 'Oligodendrocyte-Immature', 'Astrocyte-Fibrous', 'Cajal-Retzius cell', 'Oligodendrocyte', 'IPC-Glia']
        Length: 34
        Categories (34, object): ['Astrocyte-Fibrous', 'Astrocyte-Immature', 'Astrocyte-Protoplasmic',
                                'Cajal-Retzius cell', ..., 'RG-tRG', 'RG-vRG', 'Unknown', 'Vascular']
        Type_updated
        ['EN-IT-Immature', 'EN-Newborn', 'EN-L5-ET', 'RG-vRG', 'IN-MGE-Immature', ..., 'Oligodendrocyte-Immature', 'Astrocyte-Fibrous', 'Cajal-Retzius cell', 'Oligodendrocyte', 'Tri-IPC']
        Length: 34
        Categories (34, object): ['Astrocyte-Fibrous', 'Astrocyte-Immature', 'Astrocyte-Protoplasmic',
                                'Cajal-Retzius cell', ..., 'RG-vRG', 'Tri-IPC', 'Unknown', 'Vascular']
        Cluster
        ['0', '10', '78', '24', '30', ..., '91', '159', '170', '124', '96']
        Length: 176
        Categories (176, object): ['0', '1', '2', '3', ..., '172', '173', '174', '175']
        organism_ontology_term_id
        ['NCBITaxon:9606']
        Categories (1, object): ['NCBITaxon:9606']
        tissue_ontology_term_id
        ['UBERON:0000451', 'UBERON:0000411', 'UBERON:0001950', 'UBERON:0001893', 'UBERON:0001890', 'UBERON:0013541', 'UBERON:8440010', 'UBERON:0013540']
        Categories (8, object): ['UBERON:0000411', 'UBERON:0000451', 'UBERON:0001890', 'UBERON:0001893',
                                'UBERON:0013540', 'UBERON:0013541', 'UBERON:8440010', 'UBERON:0001950']
        tissue_type
        ['tissue']
        Categories (1, object): ['tissue']
        assay_ontology_term_id
        ['EFO:0030059']
        Categories (1, object): ['EFO:0030059']
        disease_ontology_term_id
        ['PATO:0000461']
        Categories (1, object): ['PATO:0000461']
        cell_type_ontology_term_id
        ['CL:4023008', 'CL:0000679', 'CL:4023041', 'CL:0013000', 'CL:4023069', ..., 'CL:4023011', 'CL:4023072', 'CL:4023059', 'CL:0000695', 'CL:0000128']
        Length: 29
        Categories (29, object): ['CL:0000128', 'CL:0000129', 'CL:0000617', 'CL:0000679', ..., 'CL:4030064',
                                'CL:4030065', 'CL:0000127', 'unknown']
        self_reported_ethnicity_ontology_term_id
        ['unknown', 'HANCESTRO:0462', 'HANCESTRO:0568', 'HANCESTRO:0005']
        Categories (4, object): ['HANCESTRO:0005', 'HANCESTRO:0462', 'HANCESTRO:0568', 'unknown']
        development_stage_ontology_term_id
        ['HsapDv:0000051', 'HsapDv:0000055', 'HsapDv:0000049', 'HsapDv:0000050', 'HsapDv:0000053', ..., 'HsapDv:0000176', 'HsapDv:0000180', 'HsapDv:0000069', 'HsapDv:0000070', 'HsapDv:0000071']
        Length: 20
        Categories (20, object): ['HsapDv:0000028', 'HsapDv:0000030', 'HsapDv:0000046', 'HsapDv:0000048', ...,
                                'HsapDv:0000176', 'HsapDv:0000178', 'HsapDv:0000180', 'HsapDv:0000182']
        sex_ontology_term_id
        ['PATO:0000384', 'PATO:0000383']
        Categories (2, object): ['PATO:0000383', 'PATO:0000384']
        donor_id
        ['ARK14', 'ARK1', 'ARK15', 'ARK25', 'ARK28', ..., 'NIH-5376', 'NIH-5554', 'NIH-5900', 'NIH-M1154', 'NIH-M2837']
        Length: 27
        Categories (27, object): ['ARK1', 'ARK5', 'ARK14', 'ARK15', ..., 'NIH-5554', 'NIH-5900',
                                'NIH-M1154', 'NIH-M2837']
        suspension_type
        ['nucleus']
        Categories (1, object): ['nucleus']
        is_primary_data
        [ True]
        cell_type
        ['intratelencephalic-projecting glutamatergic c..., 'glutamatergic neuron', 'L5 extratelencephalic projecting glutamatergi..., 'forebrain radial glial cell', 'medial ganglionic eminence derived GABAergic ..., ..., 'lamp5 GABAergic cortical interneuron', 'brain vascular cell', 'differentiation-committed oligodendrocyte pre..., 'Cajal-Retzius cell', 'oligodendrocyte']
        Length: 29
        Categories (29, object): ['oligodendrocyte', 'microglial cell', 'GABAergic neuron',
                                'glutamatergic neuron', ..., 'L5 intratelencephalic projecting glutamatergi...,
                                'L6 intratelencephalic projecting glutamatergi..., 'astrocyte', 'unknown']
        assay
        ['10x multiome']
        Categories (1, object): ['10x multiome']
        disease
        ['normal']
        Categories (1, object): ['normal']
        organism
        ['Homo sapiens']
        Categories (1, object): ['Homo sapiens']
        sex
        ['male', 'female']
        Categories (2, object): ['female', 'male']
        tissue
        ['prefrontal cortex', 'visual cortex', 'neocortex', 'telencephalon', 'forebrain', 'Brodmann (1909) area 10', 'Brodmann (1909) area 17', 'Brodmann (1909) area 9']
        Categories (8, object): ['visual cortex', 'prefrontal cortex', 'forebrain', 'telencephalon',
                                'Brodmann (1909) area 9', 'Brodmann (1909) area 10', 'Brodmann (1909) area 17',
                                'neocortex']
        self_reported_ethnicity
        ['unknown', 'British', 'African American', 'European']
        Categories (4, object): ['European', 'British', 'African American', 'unknown']
        development_stage
        ['14th week post-fertilization stage', '18th week post-fertilization stage', '12th week post-fertilization stage', '13th week post-fertilization stage', '16th week post-fertilization stage', ..., '3-month-old stage', '7-month-old stage', '32nd week post-fertilization stage', '33rd week post-fertilization stage', '34th week post-fertilization stage']
        Length: 20
        Categories (20, object): ['Carnegie stage 21', 'Carnegie stage 23',
                                '9th week post-fertilization stage', '11th week post-fertilization stage', ...,
                                '3-month-old stage', '5-month-old stage', '7-month-old stage',
                                '9-month-old stage']
        observation_joinid
        ['7%vVlx?M1t' '-sGG5QKr%@' '-gg2D3s#iP' ... 'voEXM<)>xP' '0ubK66siJb'
        'MkFy{eOCJu']
    ```
- Find developmental age and tissue combo
        - 'First_trimester'
        - 'Second_trimester'
        - 'Third_trimester'
        - 'Infancy'
        - 'Adolescence'
```
>>> age = ['First_trimester', 'Second_trimester', 'Third_trimester', 'Infancy', 'Adolescence']
>>> for i in age:
...     print(i)
...     metadata[(metadata["Group"] == i)]["tissue"].unique()
...
First_trimester
['neocortex', 'telencephalon', 'forebrain']
Categories (8, object): ['visual cortex', 'prefrontal cortex', 'forebrain', 'telencephalon',
                         'Brodmann (1909) area 9', 'Brodmann (1909) area 10', 'Brodmann (1909) area 17',
                         'neocortex']
Second_trimester
['prefrontal cortex', 'visual cortex', 'neocortex']
Categories (8, object): ['visual cortex', 'prefrontal cortex', 'forebrain', 'telencephalon',
                         'Brodmann (1909) area 9', 'Brodmann (1909) area 10', 'Brodmann (1909) area 17',
                         'neocortex']
Third_trimester
['prefrontal cortex', 'Brodmann (1909) area 10', 'Brodmann (1909) area 17']
Categories (8, object): ['visual cortex', 'prefrontal cortex', 'forebrain', 'telencephalon',
                         'Brodmann (1909) area 9', 'Brodmann (1909) area 10', 'Brodmann (1909) area 17',
                         'neocortex']
Infancy
['Brodmann (1909) area 10', 'Brodmann (1909) area 17', 'Brodmann (1909) area 9']
Categories (8, object): ['visual cortex', 'prefrontal cortex', 'forebrain', 'telencephalon',
                         'Brodmann (1909) area 9', 'Brodmann (1909) area 10', 'Brodmann (1909) area 17',
                         'neocortex']
Adolescence
['Brodmann (1909) area 17', 'Brodmann (1909) area 9']
Categories (8, object): ['visual cortex', 'prefrontal cortex', 'forebrain', 'telencephalon',
                         'Brodmann (1909) area 9', 'Brodmann (1909) area 10', 'Brodmann (1909) area 17',
                         'neocortex']
```

- In `adata.var`, gene symbol is under column `feature_name`
```
adata.var.head()
                  gene_name  feature_is_filtered       feature_name  ... feature_biotype feature_length    feature_type
ENSG00000243485         NaN                 True        MIR1302-2HG  ...            gene            623          lncRNA
ENSG00000237613         NaN                 True            FAM138A  ...            gene            888          lncRNA
ENSG00000186092         NaN                 True              OR4F5  ...            gene           2618  protein_coding
ENSG00000238009  AL627309.1                False  ENSG00000238009.6  ...            gene            629          lncRNA
ENSG00000239945  AL627309.3                False  ENSG00000239945.1  ...            gene           1319          lncRNA

[5 rows x 7 columns]

adata.var.shape
(36406, 7)
```
- There are 36,406 genes

- Count data is store under the raw.X layer

- The cell ID is a row name in `adata.obs`
```
adata.obs.head()
                                            Sample_ID  ...  observation_joinid
ARKFrozen-18-PFC_AAACAGCCAAGGTGCA-1  ARKFrozen-18-PFC  ...          7%vVlx?M1t
ARKFrozen-18-PFC_AAACAGCCAATTGAAG-1  ARKFrozen-18-PFC  ...          -sGG5QKr%@
ARKFrozen-18-PFC_AAACAGCCACTAAATC-1  ARKFrozen-18-PFC  ...          -gg2D3s#iP
ARKFrozen-18-PFC_AAACAGCCAGGTTCAC-1  ARKFrozen-18-PFC  ...          Ej2T0jem@!
ARKFrozen-18-PFC_AAACAGCCATCCATCT-1  ARKFrozen-18-PFC  ...          #@3F+6bvd(

adata.obs.shape
(232328, 39)
```
- There are 232,328 cells