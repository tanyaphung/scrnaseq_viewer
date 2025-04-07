## General information
- Link: https://cellxgene.cziscience.com/collections/f406a653-c079-4bf9-aab6-85846c27571d
- Raw counts downloaded on 2025-02-12
```
wget https://datasets.cellxgene.cziscience.com/ecd5537e-9561-43f2-9bb0-9d0b911084fe.h5ad
```

## Data exploration
### Viewing the columns of the obs.
Index(['organism_ontology_term_id', 'tissue_ontology_term_id', 'tissue_type',
       'assay_ontology_term_id', 'disease_ontology_term_id',
       'development_stage_ontology_term_id', 'sex_ontology_term_id',
       'donor_id', 'suspension_type',
       'self_reported_ethnicity_ontology_term_id',
       'cell_type_ontology_term_id', 'is_primary_data', 'cell_type', 'assay',
       'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity',
       'development_stage', 'observation_joinid'],
      dtype='object')
### View the unique values.
Find unique values for column:  organism_ontology_term_id
['NCBITaxon:9606']
Categories (1, object): ['NCBITaxon:9606']
Find unique values for column:  tissue_ontology_term_id
['UBERON:0009835', 'UBERON:0009834', 'UBERON:0001873', 'UBERON:0002421']
Categories (4, object): ['UBERON:0001873', 'UBERON:0002421', 'UBERON:0009834', 'UBERON:0009835']
Find unique values for column:  tissue_type
['tissue']
Categories (1, object): ['tissue']
Find unique values for column:  assay_ontology_term_id
['EFO:0030059']
Categories (1, object): ['EFO:0030059']
Find unique values for column:  disease_ontology_term_id
['PATO:0000461']
Categories (1, object): ['PATO:0000461']
Find unique values for column:  development_stage_ontology_term_id
['HsapDv:0000261', 'HsapDv:0000098', 'HsapDv:0000108', 'HsapDv:0000133', 'HsapDv:0000156', 'HsapDv:0000114', 'HsapDv:0000155', 'HsapDv:0000100']
Categories (8, object): ['HsapDv:0000098', 'HsapDv:0000100', 'HsapDv:0000108', 'HsapDv:0000114',
                         'HsapDv:0000133', 'HsapDv:0000155', 'HsapDv:0000156', 'HsapDv:0000261']
Find unique values for column:  sex_ontology_term_id
['PATO:0000383', 'PATO:0000384']
Categories (2, object): ['PATO:0000383', 'PATO:0000384']
Find unique values for column:  donor_id
['postnatB_3', 'postnatB_10', 'postnatB_1', 'postnatB_5', 'postnatB_2', 'postnatB_7', 'postnatB_8', 'postnatB_4', 'postnatB_6', 'postnatB_9']
Categories (10, object): ['postnatB_1', 'postnatB_2', 'postnatB_3', 'postnatB_4', ..., 'postnatB_7',
                          'postnatB_8', 'postnatB_9', 'postnatB_10']
Find unique values for column:  suspension_type
['nucleus']
Categories (1, object): ['nucleus']
Find unique values for column:  self_reported_ethnicity_ontology_term_id
['unknown']
Categories (1, object): ['unknown']
Find unique values for column:  cell_type_ontology_term_id
['CL:0000127', 'CL:0000617', 'CL:0000129', 'CL:4033054', 'CL:0000679', ..., 'CL:0000669', 'CL:0000115', 'CL:0002453', 'CL:0000128', 'CL:0000065']
Length: 11
Categories (11, object): ['CL:0000065', 'CL:0000115', 'CL:0000127', 'CL:0000128', ..., 'CL:0000679',
                          'CL:0002453', 'CL:4023051', 'CL:4033054']
Find unique values for column:  is_primary_data
[ True]
Find unique values for column:  cell_type
['astrocyte', 'GABAergic neuron', 'microglial cell', 'perivascular cell', 'glutamatergic neuron', ..., 'pericyte', 'endothelial cell', 'oligodendrocyte precursor cell', 'oligodendrocyte', 'ependymal cell']
Length: 11
Categories (11, object): ['ependymal cell', 'endothelial cell', 'astrocyte', 'oligodendrocyte', ...,
                          'glutamatergic neuron', 'oligodendrocyte precursor cell',
                          'vascular leptomeningeal cell', 'perivascular cell']
Find unique values for column:  assay
['10x multiome']
Categories (1, object): ['10x multiome']
Find unique values for column:  disease
['normal']
Categories (1, object): ['normal']
Find unique values for column:  organism
['Homo sapiens']
Categories (1, object): ['Homo sapiens']
Find unique values for column:  sex
['female', 'male']
Categories (2, object): ['female', 'male']
Find unique values for column:  tissue
['anterior cingulate cortex', 'dorsolateral prefrontal cortex', 'caudate nucleus', 'hippocampal formation']
Categories (4, object): ['caudate nucleus', 'hippocampal formation',
                         'dorsolateral prefrontal cortex', 'anterior cingulate cortex']
Find unique values for column:  self_reported_ethnicity
['unknown']
Categories (1, object): ['unknown']
Find unique values for column:  development_stage
['infant stage', '4-year-old stage', '14-year-old stage', '39-year-old stage', '62-year-old stage', '20-year-old stage', '61-year-old stage', '6-year-old stage']
Categories (8, object): ['4-year-old stage', '6-year-old stage', '14-year-old stage',
                         '20-year-old stage', '39-year-old stage', '61-year-old stage',
                         '62-year-old stage', 'infant stage']
Find unique values for column:  observation_joinid
['+JNn_x$+ro' 'c@ouCX2_EI' 'k|%PC(@G+U' ... 'd~z6XdxDr5' 'wd;{?QpY|+'
 '3I=Od0yt~#']
### View the adata.var:
                 feature_is_filtered       feature_name feature_reference  \
ENSG00000243485                False        MIR1302-2HG    NCBITaxon:9606   
ENSG00000237613                False            FAM138A    NCBITaxon:9606   
ENSG00000186092                False              OR4F5    NCBITaxon:9606   
ENSG00000238009                False  ENSG00000238009.6    NCBITaxon:9606   
ENSG00000239945                False  ENSG00000239945.1    NCBITaxon:9606   

                feature_biotype feature_length    feature_type  
ENSG00000243485            gene            623          lncRNA  
ENSG00000237613            gene            888          lncRNA  
ENSG00000186092            gene           2618  protein_coding  
ENSG00000238009            gene            629          lncRNA  
ENSG00000239945            gene           1319          lncRNA  
Size of adata.var:
(36379, 6)
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
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         1.3258771 ]
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
  0.         0.         0.         0.49240932 0.         0.
  0.         0.         0.         0.         0.         0.49240932]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.2949455 ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         1.1194378  0.         0.         0.         0.
  0.         0.         0.         0.         0.         1.1194378 ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.28671038]
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
  0.         0.         0.         0.22368486 0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.9618081  0.         0.
  0.5923382  0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.42699927]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.24750665]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.56672513 0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.69864213 0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         1.2875737 ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.37250474 0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.37250474 0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]]
Trying to print adata.X[1:25, 1:25]
  (3, 23)	1.3258771
  (6, 15)	0.49240932
  (6, 23)	0.49240932
  (7, 23)	0.2949455
  (8, 13)	1.1194378
  (8, 23)	1.1194378
  (9, 23)	0.28671038
  (12, 15)	0.22368486
  (13, 15)	0.9618081
  (13, 18)	0.5923382
  (15, 23)	0.42699927
  (16, 23)	0.24750665
  (17, 13)	0.56672513
  (19, 15)	0.69864213
  (20, 23)	1.2875737
  (23, 2)	0.37250474
  (23, 13)	0.37250474
Trying to print adata.raw.X.A[1:25, 1:25]
[[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 1.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 1. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]]
### View the adata.obs:
                            organism_ontology_term_id tissue_ontology_term_id  \
ACC_4413_AAACCGAAGTCGCAAT-1            NCBITaxon:9606          UBERON:0009835   
ACC_4413_AAACGGATCCCTCAGT-1            NCBITaxon:9606          UBERON:0009835   
ACC_4413_AAACGGATCCCTTGCG-1            NCBITaxon:9606          UBERON:0009835   
ACC_4413_AAACGTACATCGCTCC-1            NCBITaxon:9606          UBERON:0009835   
ACC_4413_AAAGCAAGTATTGGAT-1            NCBITaxon:9606          UBERON:0009835   

                            tissue_type assay_ontology_term_id  \
ACC_4413_AAACCGAAGTCGCAAT-1      tissue            EFO:0030059   
ACC_4413_AAACGGATCCCTCAGT-1      tissue            EFO:0030059   
ACC_4413_AAACGGATCCCTTGCG-1      tissue            EFO:0030059   
ACC_4413_AAACGTACATCGCTCC-1      tissue            EFO:0030059   
ACC_4413_AAAGCAAGTATTGGAT-1      tissue            EFO:0030059   

                            disease_ontology_term_id  \
ACC_4413_AAACCGAAGTCGCAAT-1             PATO:0000461   
ACC_4413_AAACGGATCCCTCAGT-1             PATO:0000461   
ACC_4413_AAACGGATCCCTTGCG-1             PATO:0000461   
ACC_4413_AAACGTACATCGCTCC-1             PATO:0000461   
ACC_4413_AAAGCAAGTATTGGAT-1             PATO:0000461   

                            development_stage_ontology_term_id  \
ACC_4413_AAACCGAAGTCGCAAT-1                     HsapDv:0000261   
ACC_4413_AAACGGATCCCTCAGT-1                     HsapDv:0000261   
ACC_4413_AAACGGATCCCTTGCG-1                     HsapDv:0000261   
ACC_4413_AAACGTACATCGCTCC-1                     HsapDv:0000261   
ACC_4413_AAAGCAAGTATTGGAT-1                     HsapDv:0000261   

                            sex_ontology_term_id    donor_id suspension_type  \
ACC_4413_AAACCGAAGTCGCAAT-1         PATO:0000383  postnatB_3         nucleus   
ACC_4413_AAACGGATCCCTCAGT-1         PATO:0000383  postnatB_3         nucleus   
ACC_4413_AAACGGATCCCTTGCG-1         PATO:0000383  postnatB_3         nucleus   
ACC_4413_AAACGTACATCGCTCC-1         PATO:0000383  postnatB_3         nucleus   
ACC_4413_AAAGCAAGTATTGGAT-1         PATO:0000383  postnatB_3         nucleus   

                            self_reported_ethnicity_ontology_term_id  \
ACC_4413_AAACCGAAGTCGCAAT-1                                  unknown   
ACC_4413_AAACGGATCCCTCAGT-1                                  unknown   
ACC_4413_AAACGGATCCCTTGCG-1                                  unknown   
ACC_4413_AAACGTACATCGCTCC-1                                  unknown   
ACC_4413_AAAGCAAGTATTGGAT-1                                  unknown   

                            cell_type_ontology_term_id  is_primary_data  \
ACC_4413_AAACCGAAGTCGCAAT-1                 CL:0000127             True   
ACC_4413_AAACGGATCCCTCAGT-1                 CL:0000617             True   
ACC_4413_AAACGGATCCCTTGCG-1                 CL:0000129             True   
ACC_4413_AAACGTACATCGCTCC-1                 CL:4033054             True   
ACC_4413_AAAGCAAGTATTGGAT-1                 CL:0000127             True   

                                     cell_type         assay disease  \
ACC_4413_AAACCGAAGTCGCAAT-1          astrocyte  10x multiome  normal   
ACC_4413_AAACGGATCCCTCAGT-1   GABAergic neuron  10x multiome  normal   
ACC_4413_AAACGGATCCCTTGCG-1    microglial cell  10x multiome  normal   
ACC_4413_AAACGTACATCGCTCC-1  perivascular cell  10x multiome  normal   
ACC_4413_AAAGCAAGTATTGGAT-1          astrocyte  10x multiome  normal   

                                 organism     sex                     tissue  \
ACC_4413_AAACCGAAGTCGCAAT-1  Homo sapiens  female  anterior cingulate cortex   
ACC_4413_AAACGGATCCCTCAGT-1  Homo sapiens  female  anterior cingulate cortex   
ACC_4413_AAACGGATCCCTTGCG-1  Homo sapiens  female  anterior cingulate cortex   
ACC_4413_AAACGTACATCGCTCC-1  Homo sapiens  female  anterior cingulate cortex   
ACC_4413_AAAGCAAGTATTGGAT-1  Homo sapiens  female  anterior cingulate cortex   

                            self_reported_ethnicity development_stage  \
ACC_4413_AAACCGAAGTCGCAAT-1                 unknown      infant stage   
ACC_4413_AAACGGATCCCTCAGT-1                 unknown      infant stage   
ACC_4413_AAACGGATCCCTTGCG-1                 unknown      infant stage   
ACC_4413_AAACGTACATCGCTCC-1                 unknown      infant stage   
ACC_4413_AAAGCAAGTATTGGAT-1                 unknown      infant stage   

                            observation_joinid  
ACC_4413_AAACCGAAGTCGCAAT-1         +JNn_x$+ro  
ACC_4413_AAACGGATCCCTCAGT-1         c@ouCX2_EI  
ACC_4413_AAACGGATCCCTTGCG-1         k|%PC(@G+U  
ACC_4413_AAACGTACATCGCTCC-1         jpq**sAdW_  
ACC_4413_AAAGCAAGTATTGGAT-1         X#t*$UV9ie  
Size of adata.obs:
(101924, 21)

## Conclusions on filtering
- by tissue: 
    - caudate nucleus
    - hippocampal formation
    - dorsolateral prefrontal cortex
    - anterior cingulate cortex
- by developmental age
    - PostnatalEarly: 
        - 'infant stage'
        - '4-year-old stage'
        - '6-year-old stage'
    - PostnatalLate: 
        - '14-year-old stage'
        - '20-year-old stage'
        - '39-year-old stage'
        - '61-year-old stage'
        - '62-year-old stage'

