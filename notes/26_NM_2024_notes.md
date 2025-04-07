## General information
- Paper: https://www.nature.com/articles/s41597-024-04117-y
- Raw counts downloaded on 2025-02-12

```
wget https://datasets.cellxgene.cziscience.com/5c56796c-849f-4919-a680-15f1a3adcd91.h5ad
```

## Data exploration

### Viewing the columns of the obs.
Index(['n_genes', 'n_counts', 'Brain_bank', 'RIN', 'path_braak_lb',
       'derived_class2', 'PMI', 'organism_ontology_term_id',
       'tissue_ontology_term_id', 'tissue_type', 'assay_ontology_term_id',
       'disease_ontology_term_id', 'cell_type_ontology_term_id',
       'self_reported_ethnicity_ontology_term_id',
       'development_stage_ontology_term_id', 'sex_ontology_term_id',
       'donor_id', 'suspension_type', 'is_primary_data', 'cell_type', 'assay',
       'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity',
       'development_stage', 'observation_joinid'],
      dtype='object')
### View the unique values.
Find unique values for column:  n_genes
[ 3633  3364  2221 ... 12040 10622 11467]
Find unique values for column:  n_counts
[10185.  9755.  4034. ... 35792. 38624. 27873.]
Find unique values for column:  Brain_bank
['MSSM', 'UM', 'Harvard']
Categories (3, object): ['Harvard', 'MSSM', 'UM']
Find unique values for column:  RIN
[4.5 4.3 3.8 2.5 3.6 4.2 3.7 4.4 4.  4.1 3.  2.6 5.1 5.6 5.  4.7 1.3 2.2
 1.2 4.8 6.3 1.4 2.8]
Find unique values for column:  path_braak_lb
[0 1 2]
Find unique values for column:  derived_class2
['OPC', 'Oligo', 'Endo', 'Micro_PVM', 'Astro', ..., 'DMNX_Neu', 'Adaptive', 'EN', 'IN', 'GPI_Neu']
Length: 12
Categories (12, object): ['DMNX_Neu' < 'GPI_Neu' < 'EN' < 'IN' ... 'Adaptive' < 'Mural' < 'Endo' <
                          'Ependymal']
Find unique values for column:  PMI
[20.91       22.         21.         15.38333333  4.55        8.5
  5.25       26.         14.          9.5        19.         17.66
         nan  8.75       17.66666667 28.         19.5         3.83333333
 22.61       30.         15.98       20.15       14.25       22.22
 24.         28.25        7.55        4.          7.         14.4
 15.         23.         19.29        4.5         4.3        18.
  8.41666667  9.7         2.91666667 20.         13.66        5.5
  7.25       25.          5.4        19.16666667 16.25       16.6
  5.33333333 27.25       30.5        18.08       11.25       15.25
 10.15       21.12       22.1        16.82       20.67       13.92
 22.88       10.6        23.5        18.83        9.8        35.42
 22.58       19.9         9.67       22.66666667 10.16666667]
Find unique values for column:  organism_ontology_term_id
['NCBITaxon:9606']
Categories (1, object): ['NCBITaxon:9606']
Find unique values for column:  tissue_ontology_term_id
['UBERON:0002870', 'UBERON:0002436', 'UBERON:0000451', 'UBERON:0001384', 'UBERON:0002477']
Categories (5, object): ['UBERON:0002870' < 'UBERON:0002477' < 'UBERON:0001384' < 'UBERON:0000451' <
                         'UBERON:0002436']
Find unique values for column:  tissue_type
['tissue']
Categories (1, object): ['tissue']
Find unique values for column:  assay_ontology_term_id
['EFO:0009922']
Categories (1, object): ['EFO:0009922']
Find unique values for column:  disease_ontology_term_id
['PATO:0000461']
Categories (1, object): ['PATO:0000461']
Find unique values for column:  cell_type_ontology_term_id
['CL:0002453', 'CL:0000128', 'CL:0000115', 'CL:0000878', 'CL:0000127', ..., 'CL:0000065', 'CL:2000029', 'CL:0000738', 'CL:0000679', 'CL:0000617']
Length: 11
Categories (11, object): ['CL:0000065', 'CL:0000115', 'CL:0000127', 'CL:0000128', ..., 'CL:0000878',
                          'CL:0002453', 'CL:0008034', 'CL:2000029']
Find unique values for column:  self_reported_ethnicity_ontology_term_id
['HANCESTRO:0005', 'HANCESTRO:0006', 'unknown', 'HANCESTRO:0010']
Categories (4, object): ['HANCESTRO:0005', 'HANCESTRO:0006', 'HANCESTRO:0010', 'unknown']
Find unique values for column:  development_stage_ontology_term_id
['HsapDv:0000162', 'HsapDv:0000170', 'HsapDv:0000166', 'HsapDv:0000160', 'HsapDv:0000168', ..., 'HsapDv:0000212', 'HsapDv:0000172', 'HsapDv:0000159', 'HsapDv:0000169', 'HsapDv:0000173']
Length: 16
Categories (16, object): ['HsapDv:0000158', 'HsapDv:0000159', 'HsapDv:0000160', 'HsapDv:0000161', ...,
                          'HsapDv:0000210', 'HsapDv:0000212', 'HsapDv:0000221', 'HsapDv:0000229']
Find unique values for column:  sex_ontology_term_id
['PATO:0000384', 'PATO:0000383']
Categories (2, object): ['PATO:0000383', 'PATO:0000384']
Find unique values for column:  donor_id
['CTRL_1', 'CTRL_2', 'CTRL_3', 'CTRL_4', 'CTRL_5', ..., 'CTRL_20', 'CTRL_21', 'CTRL_22', 'CTRL_23', 'CTRL_24']
Length: 24
Categories (24, object): ['CTRL_18', 'CTRL_19', 'CTRL_20', 'CTRL_17', ..., 'CTRL_7', 'CTRL_11',
                          'CTRL_13', 'CTRL_9']
Find unique values for column:  suspension_type
['nucleus']
Categories (1, object): ['nucleus']
Find unique values for column:  is_primary_data
[ True]
Find unique values for column:  cell_type
['oligodendrocyte precursor cell', 'oligodendrocyte', 'endothelial cell', 'central nervous system macrophage', 'astrocyte', ..., 'ependymal cell', 'central nervous system neuron', 'leukocyte', 'glutamatergic neuron', 'GABAergic neuron']
Length: 11
Categories (11, object): ['ependymal cell', 'endothelial cell', 'astrocyte', 'oligodendrocyte', ...,
                          'central nervous system macrophage', 'oligodendrocyte precursor cell', 'mural cell',
                          'central nervous system neuron']
Find unique values for column:  assay
['10x 3' v3']
Categories (1, object): ['10x 3' v3']
Find unique values for column:  disease
['normal']
Categories (1, object): ['normal']
Find unique values for column:  organism
['Homo sapiens']
Categories (1, object): ['Homo sapiens']
Find unique values for column:  sex
['male', 'female']
Categories (2, object): ['female', 'male']
Find unique values for column:  tissue
['dorsal motor nucleus of vagus nerve', 'primary visual cortex', 'prefrontal cortex', 'primary motor cortex', 'medial globus pallidus']
Categories (5, object): ['dorsal motor nucleus of vagus nerve' < 'medial globus pallidus' < 'primary motor cortex' <
                         'prefrontal cortex' < 'primary visual cortex']
Find unique values for column:  self_reported_ethnicity
['European', 'South Asian', 'unknown', 'African']
Categories (4, object): ['European', 'South Asian', 'African', 'unknown']
Find unique values for column:  development_stage
['68-year-old stage', '76-year-old stage', '72-year-old stage', '66-year-old stage', '74-year-old stage', ..., '86-year-old stage', '78-year-old stage', '65-year-old stage', '75-year-old stage', '79-year-old stage']
Length: 16
Categories (16, object): ['64-year-old stage', '65-year-old stage', '66-year-old stage',
                          '67-year-old stage', ..., '84-year-old stage', '86-year-old stage',
                          '95-year-old stage', '101-year-old stage']
Find unique values for column:  observation_joinid
['C5><2Y;jHI' 'F0`b^-><jQ' 'tu7hysS1w^' ... '71<UMtZrMw' '25FK4DT!Rl'
 '@vUY6<XF#a']
### View the adata.var:
                gene_name  n_cells  feature_is_filtered feature_name  \
gene_id
ENSG00000186827   TNFRSF4     7846                False      TNFRSF4
ENSG00000186891  TNFRSF18     9000                False     TNFRSF18
ENSG00000160072    ATAD3B   340777                False       ATAD3B
ENSG00000041988     THAP3   209923                False        THAP3
ENSG00000142611    PRDM16   283396                False       PRDM16

                feature_reference feature_biotype feature_length  \
gene_id
ENSG00000186827    NCBITaxon:9606            gene           1039
ENSG00000186891    NCBITaxon:9606            gene            789
ENSG00000160072    NCBITaxon:9606            gene           3300
ENSG00000041988    NCBITaxon:9606            gene            931
ENSG00000142611    NCBITaxon:9606            gene           3730

                   feature_type
gene_id
ENSG00000186827  protein_coding
ENSG00000186891  protein_coding
ENSG00000160072  protein_coding
ENSG00000041988  protein_coding
ENSG00000142611  protein_coding
Size of adata.var:
(17267, 8)
### Check which layer stores the raw count:
Trying to print adata.X.A[1:25, 1:25]
command adata.X.A[1:25, 1:25] failed
Trying to print adata.X[1:25, 1:25]
<Compressed Sparse Row sparse matrix of dtype 'float32'
        with 104 stored elements and shape (24, 24)>
  Coords        Values
  (0, 1)        2.420470714569092
  (0, 4)        3.737779378890991
  (0, 16)       2.420470714569092
  (2, 2)        2.3756496906280518
  (2, 9)        2.3756496906280518
  (2, 14)       2.3756496906280518
  (2, 15)       2.3756496906280518
  (2, 16)       2.3756496906280518
  (3, 4)        2.2490670680999756
  (3, 15)       2.2490670680999756
  (3, 22)       2.2490670680999756
  (4, 5)        3.2624335289001465
  (4, 12)       3.2624335289001465
  (4, 16)       3.2624335289001465
  (4, 21)       3.2624335289001465
  (5, 1)        2.4425206184387207
  (5, 5)        2.4425206184387207
  (5, 14)       3.09122371673584
  (5, 15)       2.4425206184387207
  (5, 16)       2.4425206184387207
  (6, 9)        3.56070613861084
  (6, 14)       3.56070613861084
  (6, 16)       4.239542007446289
  (6, 22)       3.56070613861084
  (7, 16)       3.019536018371582
  :     :
  (19, 12)      1.825019121170044
  (19, 14)      2.43412446975708
  (19, 16)      2.809929132461548
  (19, 17)      1.825019121170044
  (19, 22)      1.825019121170044
  (20, 4)       2.766446113586426
  (20, 5)       3.150726795196533
  (20, 9)       2.1342859268188477
  (20, 17)      2.1342859268188477
  (20, 23)      2.766446113586426
  (21, 5)       2.657763957977295
  (21, 13)      2.657763957977295
  (21, 16)      2.657763957977295
  (21, 17)      3.315229654312134
  (22, 12)      2.9992573261260986
  (22, 14)      1.9956352710723877
  (22, 16)      2.618398904800415
  (23, 2)       1.5550336837768555
  (23, 4)       1.5550336837768555
  (23, 12)      2.1365885734558105
  (23, 13)      1.5550336837768555
  (23, 14)      1.5550336837768555
  (23, 15)      1.5550336837768555
  (23, 16)      2.1365885734558105
  (23, 22)      1.5550336837768555
Trying to print adata.raw.X.A[1:25, 1:25]
command adata.raw.X.A[1:25, 1:25] failed
### View the adata.obs:
                           n_genes  n_counts Brain_bank  RIN  path_braak_lb  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT     3633   10185.0       MSSM  4.5              0
Set10_C1-AAACGCTGTGCCTGAC     3364    9755.0       MSSM  4.5              0
Set10_C1-AAAGAACAGAACGTGC     2221    4034.0       MSSM  4.5              0
Set10_C1-AAAGAACGTATGAGGC     4384   10248.0       MSSM  4.5              0
Set10_C1-AAAGAACTCAGTGCGC     4032   11795.0       MSSM  4.5              0

                          derived_class2    PMI organism_ontology_term_id  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT            OPC  20.91            NCBITaxon:9606
Set10_C1-AAACGCTGTGCCTGAC          Oligo  20.91            NCBITaxon:9606
Set10_C1-AAAGAACAGAACGTGC          Oligo  20.91            NCBITaxon:9606
Set10_C1-AAAGAACGTATGAGGC          Oligo  20.91            NCBITaxon:9606
Set10_C1-AAAGAACTCAGTGCGC          Oligo  20.91            NCBITaxon:9606

                          tissue_ontology_term_id tissue_type  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT          UBERON:0002870      tissue
Set10_C1-AAACGCTGTGCCTGAC          UBERON:0002870      tissue
Set10_C1-AAAGAACAGAACGTGC          UBERON:0002870      tissue
Set10_C1-AAAGAACGTATGAGGC          UBERON:0002870      tissue
Set10_C1-AAAGAACTCAGTGCGC          UBERON:0002870      tissue

                          assay_ontology_term_id disease_ontology_term_id  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT            EFO:0009922             PATO:0000461
Set10_C1-AAACGCTGTGCCTGAC            EFO:0009922             PATO:0000461
Set10_C1-AAAGAACAGAACGTGC            EFO:0009922             PATO:0000461
Set10_C1-AAAGAACGTATGAGGC            EFO:0009922             PATO:0000461
Set10_C1-AAAGAACTCAGTGCGC            EFO:0009922             PATO:0000461

                          cell_type_ontology_term_id  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT                 CL:0002453
Set10_C1-AAACGCTGTGCCTGAC                 CL:0000128
Set10_C1-AAAGAACAGAACGTGC                 CL:0000128
Set10_C1-AAAGAACGTATGAGGC                 CL:0000128
Set10_C1-AAAGAACTCAGTGCGC                 CL:0000128

                          self_reported_ethnicity_ontology_term_id  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT                           HANCESTRO:0005
Set10_C1-AAACGCTGTGCCTGAC                           HANCESTRO:0005
Set10_C1-AAAGAACAGAACGTGC                           HANCESTRO:0005
Set10_C1-AAAGAACGTATGAGGC                           HANCESTRO:0005
Set10_C1-AAAGAACTCAGTGCGC                           HANCESTRO:0005

                          development_stage_ontology_term_id  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT                     HsapDv:0000162
Set10_C1-AAACGCTGTGCCTGAC                     HsapDv:0000162
Set10_C1-AAAGAACAGAACGTGC                     HsapDv:0000162
Set10_C1-AAAGAACGTATGAGGC                     HsapDv:0000162
Set10_C1-AAAGAACTCAGTGCGC                     HsapDv:0000162

                          sex_ontology_term_id donor_id suspension_type  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT         PATO:0000384   CTRL_1         nucleus
Set10_C1-AAACGCTGTGCCTGAC         PATO:0000384   CTRL_1         nucleus
Set10_C1-AAAGAACAGAACGTGC         PATO:0000384   CTRL_1         nucleus
Set10_C1-AAAGAACGTATGAGGC         PATO:0000384   CTRL_1         nucleus
Set10_C1-AAAGAACTCAGTGCGC         PATO:0000384   CTRL_1         nucleus

                           is_primary_data                       cell_type  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT             True  oligodendrocyte precursor cell
Set10_C1-AAACGCTGTGCCTGAC             True                 oligodendrocyte
Set10_C1-AAAGAACAGAACGTGC             True                 oligodendrocyte
Set10_C1-AAAGAACGTATGAGGC             True                 oligodendrocyte
Set10_C1-AAAGAACTCAGTGCGC             True                 oligodendrocyte

                               assay disease      organism   sex  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT  10x 3' v3  normal  Homo sapiens  male
Set10_C1-AAACGCTGTGCCTGAC  10x 3' v3  normal  Homo sapiens  male
Set10_C1-AAAGAACAGAACGTGC  10x 3' v3  normal  Homo sapiens  male
Set10_C1-AAAGAACGTATGAGGC  10x 3' v3  normal  Homo sapiens  male
Set10_C1-AAAGAACTCAGTGCGC  10x 3' v3  normal  Homo sapiens  male

                                                        tissue  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT  dorsal motor nucleus of vagus nerve
Set10_C1-AAACGCTGTGCCTGAC  dorsal motor nucleus of vagus nerve
Set10_C1-AAAGAACAGAACGTGC  dorsal motor nucleus of vagus nerve
Set10_C1-AAAGAACGTATGAGGC  dorsal motor nucleus of vagus nerve
Set10_C1-AAAGAACTCAGTGCGC  dorsal motor nucleus of vagus nerve

                          self_reported_ethnicity  development_stage  \
barcodekey
Set10_C1-AAACGCTAGGGCCTCT                European  68-year-old stage
Set10_C1-AAACGCTGTGCCTGAC                European  68-year-old stage
Set10_C1-AAAGAACAGAACGTGC                European  68-year-old stage
Set10_C1-AAAGAACGTATGAGGC                European  68-year-old stage
Set10_C1-AAAGAACTCAGTGCGC                European  68-year-old stage

                          observation_joinid
barcodekey
Set10_C1-AAACGCTAGGGCCTCT         C5><2Y;jHI
Set10_C1-AAACGCTGTGCCTGAC         F0`b^-><jQ
Set10_C1-AAAGAACAGAACGTGC         tu7hysS1w^
Set10_C1-AAAGAACGTATGAGGC         XfNV+$wj4}
Set10_C1-AAAGAACTCAGTGCGC         Rg1n#R8++M
Size of adata.obs:
(445297, 28)

