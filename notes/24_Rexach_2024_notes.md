## General information
- Paper: https://www.sciencedirect.com/science/article/pii/S0092867424009103?via%3Dihub
- Raw counts downloaded on 2025-02-11

```
wget https://datasets.cellxgene.cziscience.com/2808a16d-64c5-451b-91a5-c1a2d9f270d5.h5ad
```

## Data exploration

### Viewing the columns of the obs.
Index(['organism_ontology_term_id', 'tissue_ontology_term_id',
       'assay_ontology_term_id', 'disease_ontology_term_id',
       'cell_type_ontology_term_id',
       'self_reported_ethnicity_ontology_term_id',
       'development_stage_ontology_term_id', 'sex_ontology_term_id',
       'donor_id', 'suspension_type', 'ct_subcluster', 'library_id',
       'tissue_type', 'is_primary_data', 'cell_type', 'assay', 'disease',
       'organism', 'sex', 'tissue', 'self_reported_ethnicity',
       'development_stage', 'observation_joinid'],
      dtype='object')
### View the unique values.
Find unique values for column:  organism_ontology_term_id
['NCBITaxon:9606']
Categories (1, object): ['NCBITaxon:9606']
Find unique values for column:  tissue_ontology_term_id
['UBERON:0002436', 'UBERON:0034891', 'UBERON:0013535']
Categories (3, object): ['UBERON:0002436', 'UBERON:0034891', 'UBERON:0013535']
Find unique values for column:  assay_ontology_term_id
['EFO:0009899', 'EFO:0009922']
Categories (2, object): ['EFO:0009899', 'EFO:0009922']
Find unique values for column:  disease_ontology_term_id
['MONDO:0019037', 'MONDO:0004975', 'MONDO:0008243', 'PATO:0000461']
Categories (4, object): ['MONDO:0004975', 'MONDO:0008243', 'MONDO:0019037', 'PATO:0000461']
Find unique values for column:  cell_type_ontology_term_id
['CL:0000679', 'CL:0000127', 'CL:0000128', 'CL:0002453', 'CL:0000498', 'CL:0000129', 'CL:0002139', 'CL:0000084', 'CL:0000669']
Categories (9, object): ['CL:0000127', 'CL:0002139', 'CL:0000679', 'CL:0000498', ..., 'CL:0000128',
                         'CL:0002453', 'CL:0000669', 'CL:0000084']
Find unique values for column:  self_reported_ethnicity_ontology_term_id
['HANCESTRO:0590']
Categories (1, object): ['HANCESTRO:0590']
Find unique values for column:  development_stage_ontology_term_id
['HsapDv:0000095', 'HsapDv:0000158', 'HsapDv:0000272', 'HsapDv:0000155', 'HsapDv:0000152', 'HsapDv:0000154', 'HsapDv:0000153', 'HsapDv:0000151', 'HsapDv:0000150']
Categories (9, object): ['HsapDv:0000095', 'HsapDv:0000150', 'HsapDv:0000151', 'HsapDv:0000152', ...,
                         'HsapDv:0000154', 'HsapDv:0000155', 'HsapDv:0000158', 'HsapDv:0000272']
Find unique values for column:  sex_ontology_term_id
['PATO:0000384', 'PATO:0000383']
Categories (2, object): ['PATO:0000383', 'PATO:0000384']
Find unique values for column:  donor_id
['R24_P2745', 'R24_P2593', 'R24_P2533', 'R24_P2734', 'R24_P2462', ..., 'R24_P2730', 'R24_P2538', 'R24_P2739', 'R24_P2321', 'R24_P2574']
Length: 40
Categories (40, object): ['R24_P2262', 'R24_P2287', 'R24_P2321', 'R24_P2323', ..., 'R24_P2742',
                          'R24_P2743', 'R24_P2744', 'R24_P2745']
Find unique values for column:  suspension_type
['nucleus']
Categories (1, object): ['nucleus']
Find unique values for column:  ct_subcluster
['V1-excitatory-1', 'V1-astrocyte-4', 'V1-excitatory-8', 'V1-excitatory-4', 'V1-excitatory-14', ..., 'BA4-astrocyte-15', 'BA4-microglia-1', 'BA4-oligodendrocyte-11', 'BA4-pericyte-4', 'BA4-opc-8']
Length: 182
Categories (182, object): ['BA4-astrocyte-0', 'BA4-astrocyte-1', 'BA4-astrocyte-2', 'BA4-astrocyte-3',
                           ..., 'V1-opc-0', 'V1-opc-1', 'V1-opc-2', 'V1-opc-3']
Find unique values for column:  library_id
['C1_1', 'C1_2', 'C1_3', 'C1_4', 'C1_5', ..., 'P5_4', 'P5_5', 'P5_6', 'P5_7', 'P5_8']
Length: 101
Categories (101, object): ['C1_1', 'C1_2', 'C1_3', 'C1_4', ..., 'P5_5', 'P5_6', 'P5_7', 'P5_8']
Find unique values for column:  tissue_type
['tissue']
Categories (1, object): ['tissue']
Find unique values for column:  is_primary_data
[ True]
Find unique values for column:  cell_type
['glutamatergic neuron', 'astrocyte', 'oligodendrocyte', 'oligodendrocyte precursor cell', 'inhibitory interneuron', 'microglial cell', 'endothelial cell of vascular tree', 'T cell', 'pericyte']
Categories (9, object): ['astrocyte', 'endothelial cell of vascular tree',
                         'glutamatergic neuron', 'inhibitory interneuron', ..., 'oligodendrocyte',
                         'oligodendrocyte precursor cell', 'pericyte', 'T cell']
Find unique values for column:  assay
['10x 3' v2', '10x 3' v3']
Categories (2, object): ['10x 3' v2', '10x 3' v3']
Find unique values for column:  disease
['progressive supranuclear palsy', 'Alzheimer disease', 'Pick disease', 'normal']
Categories (4, object): ['Alzheimer disease', 'Pick disease', 'progressive supranuclear palsy', 'normal']
Find unique values for column:  organism
['Homo sapiens']
Categories (1, object): ['Homo sapiens']
Find unique values for column:  sex
['male', 'female']
Categories (2, object): ['female', 'male']
Find unique values for column:  tissue
['primary visual cortex', 'insular cortex', 'Brodmann (1909) area 4']
Categories (3, object): ['primary visual cortex', 'insular cortex', 'Brodmann (1909) area 4']
Find unique values for column:  self_reported_ethnicity
['European American']
Categories (1, object): ['European American']
Find unique values for column:  development_stage
['80 year-old and over stage', '64-year-old stage', '60-79 year-old stage', '61-year-old stage', '58-year-old stage', '60-year-old stage', '59-year-old stage', '57-year-old stage', '56-year-old stage']
Categories (9, object): ['80 year-old and over stage', '56-year-old stage', '57-year-old stage',
                         '58-year-old stage', ..., '60-year-old stage', '61-year-old stage',
                         '64-year-old stage', '60-79 year-old stage']
Find unique values for column:  observation_joinid
['n@FliNgoSe' '?aN*YdE%Wl' 'mXSG%*e%Cx' ... '*v7Neq6Bac' 'iBWK(9A*gq'
 'e3I$L5r*@-']
### View the adata.var:
                 vst.mean  vst.variance  vst.variance.expected   
gene_ids                                                         
ENSG00000243485  0.000066      0.000066               0.000068  \
ENSG00000186092  0.000124      0.000127               0.000128   
ENSG00000238009  0.043070      0.048379               0.058816   
ENSG00000239906  0.000091      0.000091               0.000094   
ENSG00000236601  0.000011      0.000011               0.000011   

                 vst.variance.standardized  vst.variable  feature_is_filtered   
gene_ids                                                                        
ENSG00000243485                   0.972672         False                False  \
ENSG00000186092                   0.990790         False                False   
ENSG00000238009                   0.822551         False                False   
ENSG00000239906                   0.969553         False                False   
ENSG00000236601                   0.984057         False                False   

                      feature_name feature_reference feature_biotype   
gene_ids                                                               
ENSG00000243485        MIR1302-2HG    NCBITaxon:9606            gene  \
ENSG00000186092              OR4F5    NCBITaxon:9606            gene   
ENSG00000238009  ENSG00000238009.6    NCBITaxon:9606            gene   
ENSG00000239906  ENSG00000239906.1    NCBITaxon:9606            gene   
ENSG00000236601  ENSG00000236601.2    NCBITaxon:9606            gene   

                feature_length    feature_type  
gene_ids                                        
ENSG00000243485            623          lncRNA  
ENSG00000186092           2618  protein_coding  
ENSG00000238009            629          lncRNA  
ENSG00000239906            323          lncRNA  
ENSG00000236601            607          lncRNA  
Size of adata.var:
(29968, 11)
### Check which layer stores the raw count:
Trying to print adata.X.A[1:10, 1:10]
[[0.       0.       0.       0.       0.       0.       0.       0.
  0.      ]
 [0.       0.       0.       0.       0.       0.       0.       0.
  0.      ]
 [0.       0.       0.       0.       0.       0.       1.501942 0.
  0.      ]
 [0.       0.       0.       0.       0.       0.       0.       0.
  0.      ]
 [0.       0.       0.       0.       0.       0.       0.       0.
  0.      ]
 [0.       0.       0.       0.       0.       0.       0.       0.
  0.      ]
 [0.       0.       0.       0.       0.       0.       0.       0.
  0.      ]
 [0.       0.       0.       0.       0.       0.       0.       0.
  0.      ]
 [0.       0.       0.       0.       0.       0.       0.       0.
  0.      ]]
Trying to print adata.X[1:10, 1:10]
  (2, 6)	1.501942
Trying to print adata.raw.X.A[1:10, 1:10]
[[0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 1. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]]
### View the adata.obs:
                   organism_ontology_term_id tissue_ontology_term_id   
AAACCTGAGAATGTGT-1            NCBITaxon:9606          UBERON:0002436  \
AAACCTGCAATCGAAA-1            NCBITaxon:9606          UBERON:0002436   
AAACCTGCACCTGGTG-1            NCBITaxon:9606          UBERON:0002436   
AAACCTGGTAAATGAC-1            NCBITaxon:9606          UBERON:0002436   
AAACCTGGTCGTCTTC-1            NCBITaxon:9606          UBERON:0002436   

                   assay_ontology_term_id disease_ontology_term_id   
AAACCTGAGAATGTGT-1            EFO:0009899            MONDO:0019037  \
AAACCTGCAATCGAAA-1            EFO:0009899            MONDO:0019037   
AAACCTGCACCTGGTG-1            EFO:0009899            MONDO:0019037   
AAACCTGGTAAATGAC-1            EFO:0009899            MONDO:0019037   
AAACCTGGTCGTCTTC-1            EFO:0009899            MONDO:0019037   

                   cell_type_ontology_term_id   
AAACCTGAGAATGTGT-1                 CL:0000679  \
AAACCTGCAATCGAAA-1                 CL:0000127   
AAACCTGCACCTGGTG-1                 CL:0000679   
AAACCTGGTAAATGAC-1                 CL:0000679   
AAACCTGGTCGTCTTC-1                 CL:0000679   

                   self_reported_ethnicity_ontology_term_id   
AAACCTGAGAATGTGT-1                           HANCESTRO:0590  \
AAACCTGCAATCGAAA-1                           HANCESTRO:0590   
AAACCTGCACCTGGTG-1                           HANCESTRO:0590   
AAACCTGGTAAATGAC-1                           HANCESTRO:0590   
AAACCTGGTCGTCTTC-1                           HANCESTRO:0590   

                   development_stage_ontology_term_id sex_ontology_term_id   
AAACCTGAGAATGTGT-1                     HsapDv:0000095         PATO:0000384  \
AAACCTGCAATCGAAA-1                     HsapDv:0000095         PATO:0000384   
AAACCTGCACCTGGTG-1                     HsapDv:0000095         PATO:0000384   
AAACCTGGTAAATGAC-1                     HsapDv:0000095         PATO:0000384   
AAACCTGGTCGTCTTC-1                     HsapDv:0000095         PATO:0000384   

                     donor_id suspension_type    ct_subcluster library_id   
AAACCTGAGAATGTGT-1  R24_P2745         nucleus  V1-excitatory-1       C1_1  \
AAACCTGCAATCGAAA-1  R24_P2745         nucleus   V1-astrocyte-4       C1_1   
AAACCTGCACCTGGTG-1  R24_P2745         nucleus  V1-excitatory-1       C1_1   
AAACCTGGTAAATGAC-1  R24_P2745         nucleus  V1-excitatory-8       C1_1   
AAACCTGGTCGTCTTC-1  R24_P2745         nucleus  V1-excitatory-4       C1_1   

                   tissue_type  is_primary_data             cell_type   
AAACCTGAGAATGTGT-1      tissue             True  glutamatergic neuron  \
AAACCTGCAATCGAAA-1      tissue             True             astrocyte   
AAACCTGCACCTGGTG-1      tissue             True  glutamatergic neuron   
AAACCTGGTAAATGAC-1      tissue             True  glutamatergic neuron   
AAACCTGGTCGTCTTC-1      tissue             True  glutamatergic neuron   

                        assay                         disease      organism   
AAACCTGAGAATGTGT-1  10x 3' v2  progressive supranuclear palsy  Homo sapiens  \
AAACCTGCAATCGAAA-1  10x 3' v2  progressive supranuclear palsy  Homo sapiens   
AAACCTGCACCTGGTG-1  10x 3' v2  progressive supranuclear palsy  Homo sapiens   
AAACCTGGTAAATGAC-1  10x 3' v2  progressive supranuclear palsy  Homo sapiens   
AAACCTGGTCGTCTTC-1  10x 3' v2  progressive supranuclear palsy  Homo sapiens   

                     sex                 tissue self_reported_ethnicity   
AAACCTGAGAATGTGT-1  male  primary visual cortex       European American  \
AAACCTGCAATCGAAA-1  male  primary visual cortex       European American   
AAACCTGCACCTGGTG-1  male  primary visual cortex       European American   
AAACCTGGTAAATGAC-1  male  primary visual cortex       European American   
AAACCTGGTCGTCTTC-1  male  primary visual cortex       European American   

                             development_stage observation_joinid  
AAACCTGAGAATGTGT-1  80 year-old and over stage         n@FliNgoSe  
AAACCTGCAATCGAAA-1  80 year-old and over stage         ?aN*YdE%Wl  
AAACCTGCACCTGGTG-1  80 year-old and over stage         mXSG%*e%Cx  
AAACCTGGTAAATGAC-1  80 year-old and over stage         1m^mYbTK&X  
AAACCTGGTCGTCTTC-1  80 year-old and over stage         HClVP!<VWs  
Size of adata.obs:
(432555, 23)
