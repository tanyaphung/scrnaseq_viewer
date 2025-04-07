## General information
- Link: https://cellxgene.cziscience.com/collections/5d445965-6f1a-4b68-ba3a-b8f765155d3a
- Raw counts downloaded on 2025-04-01
```
wget https://datasets.cellxgene.cziscience.com/a874b4da-c6a8-486c-844c-71baa549c1e4.h5ad
```

### Viewing the columns of the obs.
Index(['nGene', 'nUMI', 'channel', 'region', 'percent.ribo', 'free_annotation',
       'donor_id', 'sample', 'location', 'magnetic.selection',
       'preparation.site', 'compartment', 'tissue_ontology_term_id',
       'assay_ontology_term_id', 'disease_ontology_term_id',
       'development_stage_ontology_term_id', 'cell_type_ontology_term_id',
       'self_reported_ethnicity_ontology_term_id', 'sex_ontology_term_id',
       'is_primary_data', 'organism_ontology_term_id', 'suspension_type',
       'tissue_type', 'cell_type', 'assay', 'disease', 'organism', 'sex',
       'tissue', 'self_reported_ethnicity', 'development_stage',
       'observation_joinid'],
      dtype='object')
### View the unique values.
Find unique values for column:  nGene
[1347 1713 1185 ... 4986 6109 7790]
Find unique values for column:  nUMI
[ 2914  4226  2152 ... 12619  7687  8272]
Find unique values for column:  channel
['P2_1', 'P2_2', 'P2_3', 'P2_4', 'P2_5', ..., 'P3_5', 'P3_6', 'P3_7', 'P3_8', 'P3_4']
Length: 19
Categories (19, object): ['P1_1', 'P1_2', 'P1_3', 'P1_4', ..., 'P3_5', 'P3_6', 'P3_7', 'P3_8']
Find unique values for column:  region
['normal']
Categories (1, object): ['normal']
Find unique values for column:  percent.ribo
[0.0353466  0.06105064 0.04832714 ... 0.17850953 0.14132187 0.07853188]
Find unique values for column:  free_annotation
['Capillary Aerocyte', 'Capillary', 'Capillary Intermediate 1', 'Capillary Intermediate 2', 'IGSF21+ Dendritic', ..., 'Fibromyocyte', 'Ionocyte', 'Serous', 'Proximal Ciliated', 'Goblet']
Length: 57
Categories (57, object): ['Adventitial Fibroblast', 'Airway Smooth Muscle', 'Alveolar Epithelial Type 1',
                          'Alveolar Epithelial Type 2', ..., 'Signaling Alveolar Epithelial Type 2',
                          'TREM2+ Dendritic', 'Vascular Smooth Muscle', 'Vein']
Find unique values for column:  donor_id
['2', '1', '3']
Categories (3, object): ['1', '2', '3']
Find unique values for column:  sample
['distal 2', 'medial 2', 'blood 1', 'distal 1a', 'proximal 3', 'distal 3', 'blood 3']
Categories (7, object): ['blood 1', 'blood 3', 'distal 1a', 'distal 2', 'distal 3', 'medial 2',
                         'proximal 3']
Find unique values for column:  location
['distal', 'medial', 'blood', 'proximal']
Categories (4, object): ['blood', 'distal', 'medial', 'proximal']
Find unique values for column:  magnetic.selection
['epithelial', 'immune and endothelial', 'stromal', 'blood']
Categories (4, object): ['blood', 'epithelial', 'immune and endothelial', 'stromal']
Find unique values for column:  preparation.site
['biohub']
Categories (1, object): ['biohub']
Find unique values for column:  compartment
['endothelial', 'immune', 'stromal', 'epithelial']
Categories (4, object): ['endothelial', 'epithelial', 'immune', 'stromal']
Find unique values for column:  tissue_ontology_term_id
['UBERON:0002048', 'UBERON:0000178']
Categories (2, object): ['UBERON:0000178', 'UBERON:0002048']
Find unique values for column:  assay_ontology_term_id
['EFO:0009899']
Categories (1, object): ['EFO:0009899']
Find unique values for column:  disease_ontology_term_id
['PATO:0000461']
Categories (1, object): ['PATO:0000461']
Find unique values for column:  development_stage_ontology_term_id
['HsapDv:0000140', 'HsapDv:0000169', 'HsapDv:0000145']
Categories (3, object): ['HsapDv:0000140', 'HsapDv:0000145', 'HsapDv:0000169']
Find unique values for column:  cell_type_ontology_term_id
['CL:0000115', 'CL:0002144', 'CL:0000451', 'CL:0001057', 'CL:0001058', ..., 'CL:1000223', 'CL:0017000', 'CL:0019001', 'CL:0000064', 'CL:1000143']
Length: 46
Categories (46, object): ['CL:0000057', 'CL:0000064', 'CL:0000066', 'CL:0000115', ..., 'CL:1000223',
                          'CL:1000271', 'CL:1000413', 'CL:1000491']
Find unique values for column:  self_reported_ethnicity_ontology_term_id
['unknown']
Categories (1, object): ['unknown']
Find unique values for column:  sex_ontology_term_id
['PATO:0000384', 'PATO:0000383']
Categories (2, object): ['PATO:0000383', 'PATO:0000384']
Find unique values for column:  is_primary_data
[ True]
Find unique values for column:  organism_ontology_term_id
['NCBITaxon:9606']
Categories (1, object): ['NCBITaxon:9606']
Find unique values for column:  suspension_type
['cell']
Categories (1, object): ['cell']
Find unique values for column:  tissue_type
['tissue']
Categories (1, object): ['tissue']
Find unique values for column:  cell_type
['endothelial cell', 'capillary endothelial cell', 'dendritic cell', 'myeloid dendritic cell, human', 'plasmacytoid dendritic cell, human', ..., 'lung neuroendocrine cell', 'pulmonary ionocyte', 'tracheobronchial serous cell', 'ciliated cell', 'lung goblet cell']
Length: 46
Categories (46, object): ['fibroblast', 'ciliated cell', 'epithelial cell', 'endothelial cell', ...,
                          'lung neuroendocrine cell', 'lung ciliated cell', 'endothelial cell of artery',
                          'mesothelial cell of pleura']
Find unique values for column:  assay
['10x 3' v2']
Categories (1, object): ['10x 3' v2']
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
['lung', 'blood']
Categories (2, object): ['blood', 'lung']
Find unique values for column:  self_reported_ethnicity
['unknown']
Categories (1, object): ['unknown']
Find unique values for column:  development_stage
['46-year-old stage', '75-year-old stage', '51-year-old stage']
Categories (3, object): ['46-year-old stage', '51-year-old stage', '75-year-old stage']
Find unique values for column:  observation_joinid
['`D(skAXifB' 'lkR(zTsYD^' '^(|@L(t<^&' ... 'jttJ%sZMh#' 'esp8j{jjP5'
 'pfmoC_q0&H']
### View the adata.var:
                 feature_is_filtered              feature_name  \
ENSG00000119048                False                     UBE2B   
ENSG00000162852                False                      CNST   
ENSG00000265787                False  CYP4F35P_ENSG00000265787   
ENSG00000205864                False                  KRTAP5-6   
ENSG00000185052                False                   SLC24A3   

                feature_reference feature_biotype feature_length  \
ENSG00000119048    NCBITaxon:9606            gene            669   
ENSG00000162852    NCBITaxon:9606            gene           2594   
ENSG00000265787    NCBITaxon:9606            gene            390   
ENSG00000205864    NCBITaxon:9606            gene            561   
ENSG00000185052    NCBITaxon:9606            gene           3922   

                                       feature_type  
ENSG00000119048                      protein_coding  
ENSG00000162852                      protein_coding  
ENSG00000265787  transcribed_unprocessed_pseudogene  
ENSG00000205864                      protein_coding  
ENSG00000185052                      protein_coding  
Size of adata.var:
(24835, 6)
### Check which layer stores the raw count:
Trying to print adata.X.A[1:25, 1:25]
[[0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         1.3673235  0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [1.604973   0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         1.604973   0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.96115947 0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         1.7888488  0.        ]
 [1.4196758  0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         1.4766784  0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         1.2119887
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         1.8575985  0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  1.9234506  0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  1.1194379  0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         1.1194379  0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         1.3838991  0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         1.9086273
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         2.1239092
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         1.1821628  0.
  0.         0.         0.         1.1821628  0.         0.        ]]
Trying to print adata.X[1:25, 1:25]
  (2, 14)	1.3673235
  (3, 0)	1.604973
  (3, 21)	1.604973
  (4, 6)	0.96115947
  (5, 22)	1.7888488
  (6, 0)	1.4196758
  (7, 16)	1.4766784
  (9, 11)	1.2119887
  (10, 22)	1.8575985
  (11, 12)	1.9234506
  (14, 6)	1.1194379
  (14, 21)	1.1194379
  (17, 22)	1.3838991
  (20, 11)	1.9086273
  (22, 11)	2.1239092
  (23, 16)	1.1821628
  (23, 21)	1.1821628
Trying to print adata.raw.X.A[1:25, 1:25]
[[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
 [0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 1. 0. 0.]]
### View the adata.obs:
                       nGene  nUMI channel  region  percent.ribo  \
index                                                              
P2_1_AAACCTGAGAAACCAT   1347  2914    P2_1  normal      0.035347   
P2_1_AAATGCCAGATGAGAG   1713  4226    P2_1  normal      0.061051   
P2_1_AACACGTTCGATCCCT   1185  2152    P2_1  normal      0.048327   
P2_1_AACACGTTCGCACTCT   1378  3419    P2_1  normal      0.032758   
P2_1_AACCATGCAGCTCGCA   1210  2514    P2_1  normal      0.050119   

                          free_annotation donor_id    sample location  \
index                                                                   
P2_1_AAACCTGAGAAACCAT  Capillary Aerocyte        2  distal 2   distal   
P2_1_AAATGCCAGATGAGAG  Capillary Aerocyte        2  distal 2   distal   
P2_1_AACACGTTCGATCCCT  Capillary Aerocyte        2  distal 2   distal   
P2_1_AACACGTTCGCACTCT  Capillary Aerocyte        2  distal 2   distal   
P2_1_AACCATGCAGCTCGCA  Capillary Aerocyte        2  distal 2   distal   

                      magnetic.selection preparation.site  compartment  \
index                                                                    
P2_1_AAACCTGAGAAACCAT         epithelial           biohub  endothelial   
P2_1_AAATGCCAGATGAGAG         epithelial           biohub  endothelial   
P2_1_AACACGTTCGATCCCT         epithelial           biohub  endothelial   
P2_1_AACACGTTCGCACTCT         epithelial           biohub  endothelial   
P2_1_AACCATGCAGCTCGCA         epithelial           biohub  endothelial   

                      tissue_ontology_term_id assay_ontology_term_id  \
index                                                                  
P2_1_AAACCTGAGAAACCAT          UBERON:0002048            EFO:0009899   
P2_1_AAATGCCAGATGAGAG          UBERON:0002048            EFO:0009899   
P2_1_AACACGTTCGATCCCT          UBERON:0002048            EFO:0009899   
P2_1_AACACGTTCGCACTCT          UBERON:0002048            EFO:0009899   
P2_1_AACCATGCAGCTCGCA          UBERON:0002048            EFO:0009899   

                      disease_ontology_term_id  \
index                                            
P2_1_AAACCTGAGAAACCAT             PATO:0000461   
P2_1_AAATGCCAGATGAGAG             PATO:0000461   
P2_1_AACACGTTCGATCCCT             PATO:0000461   
P2_1_AACACGTTCGCACTCT             PATO:0000461   
P2_1_AACCATGCAGCTCGCA             PATO:0000461   

                      development_stage_ontology_term_id  \
index                                                      
P2_1_AAACCTGAGAAACCAT                     HsapDv:0000140   
P2_1_AAATGCCAGATGAGAG                     HsapDv:0000140   
P2_1_AACACGTTCGATCCCT                     HsapDv:0000140   
P2_1_AACACGTTCGCACTCT                     HsapDv:0000140   
P2_1_AACCATGCAGCTCGCA                     HsapDv:0000140   

                      cell_type_ontology_term_id  \
index                                              
P2_1_AAACCTGAGAAACCAT                 CL:0000115   
P2_1_AAATGCCAGATGAGAG                 CL:0000115   
P2_1_AACACGTTCGATCCCT                 CL:0000115   
P2_1_AACACGTTCGCACTCT                 CL:0000115   
P2_1_AACCATGCAGCTCGCA                 CL:0000115   

                      self_reported_ethnicity_ontology_term_id  \
index                                                            
P2_1_AAACCTGAGAAACCAT                                  unknown   
P2_1_AAATGCCAGATGAGAG                                  unknown   
P2_1_AACACGTTCGATCCCT                                  unknown   
P2_1_AACACGTTCGCACTCT                                  unknown   
P2_1_AACCATGCAGCTCGCA                                  unknown   

                      sex_ontology_term_id  is_primary_data  \
index                                                         
P2_1_AAACCTGAGAAACCAT         PATO:0000384             True   
P2_1_AAATGCCAGATGAGAG         PATO:0000384             True   
P2_1_AACACGTTCGATCCCT         PATO:0000384             True   
P2_1_AACACGTTCGCACTCT         PATO:0000384             True   
P2_1_AACCATGCAGCTCGCA         PATO:0000384             True   

                      organism_ontology_term_id suspension_type tissue_type  \
index                                                                         
P2_1_AAACCTGAGAAACCAT            NCBITaxon:9606            cell      tissue   
P2_1_AAATGCCAGATGAGAG            NCBITaxon:9606            cell      tissue   
P2_1_AACACGTTCGATCCCT            NCBITaxon:9606            cell      tissue   
P2_1_AACACGTTCGCACTCT            NCBITaxon:9606            cell      tissue   
P2_1_AACCATGCAGCTCGCA            NCBITaxon:9606            cell      tissue   

                              cell_type      assay disease      organism  \
index                                                                      
P2_1_AAACCTGAGAAACCAT  endothelial cell  10x 3' v2  normal  Homo sapiens   
P2_1_AAATGCCAGATGAGAG  endothelial cell  10x 3' v2  normal  Homo sapiens   
P2_1_AACACGTTCGATCCCT  endothelial cell  10x 3' v2  normal  Homo sapiens   
P2_1_AACACGTTCGCACTCT  endothelial cell  10x 3' v2  normal  Homo sapiens   
P2_1_AACCATGCAGCTCGCA  endothelial cell  10x 3' v2  normal  Homo sapiens   

                        sex tissue self_reported_ethnicity  development_stage  \
index                                                                           
P2_1_AAACCTGAGAAACCAT  male   lung                 unknown  46-year-old stage   
P2_1_AAATGCCAGATGAGAG  male   lung                 unknown  46-year-old stage   
P2_1_AACACGTTCGATCCCT  male   lung                 unknown  46-year-old stage   
P2_1_AACACGTTCGCACTCT  male   lung                 unknown  46-year-old stage   
P2_1_AACCATGCAGCTCGCA  male   lung                 unknown  46-year-old stage   

                      observation_joinid  
index                                     
P2_1_AAACCTGAGAAACCAT         `D(skAXifB  
P2_1_AAATGCCAGATGAGAG         lkR(zTsYD^  
P2_1_AACACGTTCGATCCCT         ^(|@L(t<^&  
P2_1_AACACGTTCGCACTCT         !$p`uS3q4m  
P2_1_AACCATGCAGCTCGCA         z!;cZTIe`>  
Size of adata.obs:
(65662, 32)

## Data processing
1. Separate by `tissue`
- blood 
- lung
2. 1 level of cell type from the column `cell_type`
3. Note that I process only the file from the 10x 3' v2 assay (and not the Smart-seq2 assay)