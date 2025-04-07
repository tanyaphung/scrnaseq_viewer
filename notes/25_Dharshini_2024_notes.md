## General information
- Paper: https://www.biorxiv.org/content/10.1101/2024.11.03.621787v1.full.pdf
- Raw counts downloaded on 2025-02-11

```
wget https://datasets.cellxgene.cziscience.com/535e37d4-1f8f-41f2-8a3e-dd00709e6b36.h5ad
```

## Data exploration

### Viewing the columns of the obs.
Index(['APOE_class', 'Brain.Region', 'SORT', 'Braak.stage', 'Disease.Group',
       'Amyloid', 'Brain.weight', 'PMI.hr.', 'Race', 'Age', 'RIN',
       'total_counts', 'pct_counts_mt', 'droplet_FDR', 'doublefinder',
       'n_genes', 'sample', 'tissue_ontology_term_id',
       'assay_ontology_term_id', 'cell_type_ontology_term_id',
       'development_stage_ontology_term_id',
       'self_reported_ethnicity_ontology_term_id', 'disease_ontology_term_id',
       'sex_ontology_term_id', 'is_primary_data', 'organism_ontology_term_id',
       'suspension_type', 'donor_id', 'Age_group', 'PMI_group', 'RIN_group',
       'Major_celltypes', 'Author_Annotation', 'NP.Diagonis', 'tissue_type',
       'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue',
       'self_reported_ethnicity', 'development_stage', 'observation_joinid'],
      dtype='object')
### View the unique values.
Find unique values for column:  APOE_class
['E2/E4', 'E3/E3', 'E3/E4', 'E4/E4', 'nan']
Categories (5, object): ['E2/E4', 'E3/E3', 'E3/E4', 'E4/E4', 'nan']
Find unique values for column:  Brain.Region
['Frontal Cx (BA9)', 'Precuneous (BA7)', 'Primary Visual Cx (BA17)']
Categories (3, object): ['Frontal Cx (BA9)', 'Precuneous (BA7)', 'Primary Visual Cx (BA17)']
Find unique values for column:  SORT
['All Nuclei', 'NeuN']
Categories (2, object): ['All Nuclei', 'NeuN']
Find unique values for column:  Braak.stage
['IV', 'II', 'III', '0', 'VI', 'V', 'I']
Categories (7, object): ['0', 'I', 'II', 'III', 'IV', 'V', 'VI']
Find unique values for column:  Disease.Group
['int', 'low', 'high']
Categories (3, object): ['high', 'int', 'low']
Find unique values for column:  Amyloid
['DP (C0)', 'No', 'C1', 'C3', 'C2']
Categories (5, object): ['C1', 'C2', 'C3', 'DP (C0)', 'No']
Find unique values for column:  Brain.weight
['1370.0', '1480.0', '970.0', '1040.0', 'nan', ..., '980.0', '1100.0', '1300.0', '1030.0', '850.0']
Length: 27
Categories (27, object): ['850.0', '900.0', '950.0', '960.0', ..., '1400.0', '1450.0', '1480.0', 'nan']
Find unique values for column:  PMI.hr.
['unknown', '19.5', '11.8', '9.3', '10.0', ..., '10.1', '15.2', '17.0', '19.25', '10.4']
Length: 35
Categories (35, object): ['1.0', '3.5', '3.75', '4.0', ..., '24.9', '33.0', '44.0', 'unknown']
Find unique values for column:  Race
['nan', 'black', 'white', 'hispanic']
Categories (4, object): ['black', 'hispanic', 'nan', 'white']
Find unique values for column:  Age
['71', '61', '67', '87', '81', ..., '93', '82', '75', '76', '86']
Length: 27
Categories (27, object): ['52', '57', '58', '61', ..., '91', '92', '93', '95 or above']
Find unique values for column:  RIN
['4.7', '6.5', '7.8', '7.5', '6.1', ..., '5.1', '7.2', '5.4', '5.0', '6.0']
Length: 32
Categories (32, object): ['4.2', '4.3', '4.5', '4.6', ..., '7.5', '7.8', '8.6', 'unknown']
Find unique values for column:  total_counts
[4184.6753 2703.646  2393.7102 ... 3854.8196 3130.0535 2825.2693]
Find unique values for column:  pct_counts_mt
[0.14477587 0.565152   0.75240487 ... 1.154158   0.12379766 0.5450244 ]
Find unique values for column:  droplet_FDR
[0.0001358  0.         0.00012351 ... 0.04874078 0.01941958 0.00046188]
Find unique values for column:  doublefinder
['0']
Categories (1, object): ['0']
Find unique values for column:  n_genes
[3272.  999.  729. ... 5592. 6526. 5745.]
Find unique values for column:  sample
['C0060', 'C0018', 'C0023', 'C0024', 'C0025', ..., 'D0196', 'D0197', 'D0199', 'D0201', 'D0205']
Length: 243
Categories (243, object): ['C0017', 'C0018', 'C0019', 'C0020', ..., 'D0204', 'D0205', 'D0206', 'D0207']
Find unique values for column:  tissue_ontology_term_id
['UBERON:0013540', 'UBERON:0013538', 'UBERON:8440010']
Categories (3, object): ['UBERON:0013538', 'UBERON:0013540', 'UBERON:8440010']
Find unique values for column:  assay_ontology_term_id
['EFO:0009922', 'EFO:0009899', 'EFO:0008722']
Categories (3, object): ['EFO:0008722', 'EFO:0009899', 'EFO:0009922']
Find unique values for column:  cell_type_ontology_term_id
['CL:0000127', 'CL:0002453', 'CL:0000129', 'CL:0000128', 'CL:0000071', 'CL:0000679', 'CL:0000617']
Categories (7, object): ['CL:0000071', 'CL:0000127', 'CL:0000128', 'CL:0000129', 'CL:0000617',
                         'CL:0000679', 'CL:0002453']
Find unique values for column:  development_stage_ontology_term_id
['HsapDv:0000165', 'HsapDv:0000155', 'HsapDv:0000161', 'HsapDv:0000213', 'HsapDv:0000207', ..., 'HsapDv:0000252', 'HsapDv:0000208', 'HsapDv:0000169', 'HsapDv:0000170', 'HsapDv:0000212']
Length: 28
Categories (28, object): ['HsapDv:0000146', 'HsapDv:0000151', 'HsapDv:0000152', 'HsapDv:0000155', ...,
                          'HsapDv:0000218', 'HsapDv:0000219', 'HsapDv:0000223', 'HsapDv:0000252']
Find unique values for column:  self_reported_ethnicity_ontology_term_id
['unknown', 'HANCESTRO:0568', 'HANCESTRO:0005', 'HANCESTRO:0014']
Categories (4, object): ['HANCESTRO:0005', 'HANCESTRO:0014', 'HANCESTRO:0568', 'unknown']
Find unique values for column:  disease_ontology_term_id
['MONDO:0004975', 'PATO:0000461']
Categories (2, object): ['MONDO:0004975', 'PATO:0000461']
Find unique values for column:  sex_ontology_term_id
['PATO:0000384', 'PATO:0000383']
Categories (2, object): ['PATO:0000383', 'PATO:0000384']
Find unique values for column:  is_primary_data
[ True False]
Find unique values for column:  organism_ontology_term_id
['NCBITaxon:9606']
Categories (1, object): ['NCBITaxon:9606']
Find unique values for column:  suspension_type
['nucleus']
Categories (1, object): ['nucleus']
Find unique values for column:  donor_id
['19', '50', '45', '49', '30', ..., '108', '33', '95', '110', '123']
Length: 46
Categories (46, object): ['12', '14', '18', '19', ..., '137', '140', '145', '146']
Find unique values for column:  Age_group
['70 to 80', 'less than 70', '81 to 90', 'More than 90', '95 or above']
Categories (5, object): ['70 to 80', '81 to 90', '95 or above', 'More than 90', 'less than 70']
Find unique values for column:  PMI_group
['unknown', '15.2 to 19.5 hours', '10 to 15 hours', '5 to 9.3 hours', 'More than 20 hours', '1 to 4.5 hours']
Categories (6, object): ['1 to 4.5 hours', '10 to 15 hours', '15.2 to 19.5 hours', '5 to 9.3 hours',
                         'More than 20 hours', 'unknown']
Find unique values for column:  RIN_group
['4.2 to 5', '6.1 to 7', 'More than 7', '5.1 to 6', 'Unknown']
Categories (5, object): ['4.2 to 5', '5.1 to 6', '6.1 to 7', 'More than 7', 'Unknown']
Find unique values for column:  Major_celltypes
['Astrocytes', 'OPC', 'Microglia', 'Oligodendrocytes', 'Vascular', 'Excitatory', 'Inhibitory']
Categories (7, object): ['Astrocytes', 'Excitatory', 'Inhibitory', 'Microglia', 'OPC',
                         'Oligodendrocytes', 'Vascular']
Find unique values for column:  Author_Annotation
['Astrocyte-GFAP-OSMR', 'Astrocyte-SLC1A2-WIF1', 'Astrocyte-SLC1A2-SMTN', 'Astrocyte-GFAP-VCAN', 'OPC', ..., 'In5:[LHX6-SST-SPON1]', 'In8:[LHX6-SST-SEL1L3]', 'In12:[ADARB2-LAMP5-NDNF]', 'In13:[ADARB2-SYT6]', 'In10:[LHX6-ADARB2-LAMP5-HCRTR2]']
Length: 51
Categories (51, object): ['Astrocyte-GFAP-OSMR', 'Astrocyte-GFAP-VCAN', 'Astrocyte-SLC1A2-SMTN',
                          'Astrocyte-SLC1A2-WIF1', ..., 'Oligodendrocyte-COL18A1', 'Oligodendrocyte-OPALIN',
                          'Pericytes', 'VLMC']
Find unique values for column:  NP.Diagonis
['AD', 'PART', 'Control']
Categories (3, object): ['AD', 'Control', 'PART']
Find unique values for column:  tissue_type
['tissue']
Categories (1, object): ['tissue']
Find unique values for column:  cell_type
['astrocyte', 'oligodendrocyte precursor cell', 'microglial cell', 'oligodendrocyte', 'blood vessel endothelial cell', 'glutamatergic neuron', 'GABAergic neuron']
Categories (7, object): ['blood vessel endothelial cell', 'astrocyte', 'oligodendrocyte', 'microglial cell',
                         'GABAergic neuron', 'glutamatergic neuron',
                         'oligodendrocyte precursor cell']
Find unique values for column:  assay
['10x 3' v3', '10x 3' v2', 'Drop-seq']
Categories (3, object): ['Drop-seq', '10x 3' v2', '10x 3' v3']
Find unique values for column:  disease
['Alzheimer disease', 'normal']
Categories (2, object): ['Alzheimer disease', 'normal']
Find unique values for column:  organism
['Homo sapiens']
Categories (1, object): ['Homo sapiens']
Find unique values for column:  sex
['male', 'female']
Categories (2, object): ['female', 'male']
Find unique values for column:  tissue
['Brodmann (1909) area 9', 'Brodmann (1909) area 7', 'Brodmann (1909) area 17']
Categories (3, object): ['Brodmann (1909) area 7', 'Brodmann (1909) area 9', 'Brodmann (1909) area 17']
Find unique values for column:  self_reported_ethnicity
['unknown', 'African American', 'European', 'Hispanic or Latin American']
Categories (4, object): ['European', 'Hispanic or Latin American', 'African American',
                         'unknown']
Find unique values for column:  development_stage
['71-year-old stage', '61-year-old stage', '67-year-old stage', '87-year-old stage', '81-year-old stage', ..., '111-year-old stage', '82-year-old stage', '75-year-old stage', '76-year-old stage', '86-year-old stage']
Length: 28
Categories (28, object): ['52-year-old stage', '57-year-old stage', '58-year-old stage',
                          '61-year-old stage', ..., '92-year-old stage', '93-year-old stage',
                          '97-year-old stage', '111-year-old stage']
Find unique values for column:  observation_joinid
[';H|J`@{xnY' 'JQP1;<W}tH' 'Q3`kwi_dnf' ... 'MvQ=Al?cqe' '^-muzl|$EG'
 '+&08k9N-lw']
### View the adata.var:
                    mt   ribo     hb  feature_is_filtered       feature_name   
ENSEMBLE GENE                                                                  
ENSG00000186827  False  False  False                False            TNFRSF4  \
ENSG00000186891  False  False  False                False           TNFRSF18   
ENSG00000160072  False  False  False                False             ATAD3B   
ENSG00000260179  False  False  False                False  ENSG00000260179.1   
ENSG00000234396  False  False  False                False  ENSG00000234396.3   

                feature_reference feature_biotype feature_length   
ENSEMBLE GENE                                                      
ENSG00000186827    NCBITaxon:9606            gene           1039  \
ENSG00000186891    NCBITaxon:9606            gene            789   
ENSG00000160072    NCBITaxon:9606            gene           3300   
ENSG00000260179    NCBITaxon:9606            gene           1558   
ENSG00000234396    NCBITaxon:9606            gene            326   

                   feature_type  
ENSEMBLE GENE                    
ENSG00000186827  protein_coding  
ENSG00000186891  protein_coding  
ENSG00000160072  protein_coding  
ENSG00000260179          lncRNA  
ENSG00000234396          lncRNA  
Size of adata.var:
(61427, 9)
### Check which layer stores the raw count:
Trying to print adata.X.A[1:10, 1:10]
[[0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]]
Trying to print adata.X[1:10, 1:10]

Trying to print adata.raw.X.A[1:10, 1:10]
[[0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0.]]
### View the adata.obs:
                       APOE_class      Brain.Region        SORT Braak.stage   
AGTACCATCCATCTGC-C0060      E2/E4  Frontal Cx (BA9)  All Nuclei          IV  \
ACATCAGTCTTTAGTC-C0018      E3/E3  Frontal Cx (BA9)  All Nuclei          II   
CCAATCCAGAATAGGG-C0018      E3/E3  Frontal Cx (BA9)  All Nuclei          II   
CTAGCCTAGCGTTTAC-C0018      E3/E3  Frontal Cx (BA9)  All Nuclei          II   
CTGAAACGTCCCTACT-C0018      E3/E3  Frontal Cx (BA9)  All Nuclei          II   

                       Disease.Group  Amyloid Brain.weight  PMI.hr. Race Age   
AGTACCATCCATCTGC-C0060           int  DP (C0)       1370.0  unknown  nan  71  \
ACATCAGTCTTTAGTC-C0018           low       No       1480.0     19.5  nan  61   
CCAATCCAGAATAGGG-C0018           low       No       1480.0     19.5  nan  61   
CTAGCCTAGCGTTTAC-C0018           low       No       1480.0     19.5  nan  61   
CTGAAACGTCCCTACT-C0018           low       No       1480.0     19.5  nan  61   

                        RIN  total_counts  pct_counts_mt  droplet_FDR   
AGTACCATCCATCTGC-C0060  4.7     4184.6753       0.144776     0.000136  \
ACATCAGTCTTTAGTC-C0018  6.5     2703.6460       0.565152     0.000000   
CCAATCCAGAATAGGG-C0018  6.5     2393.7102       0.752405     0.000000   
CTAGCCTAGCGTTTAC-C0018  6.5     2582.0920       0.106748     0.000000   
CTGAAACGTCCCTACT-C0018  6.5     2497.0735       0.352070     0.000000   

                       doublefinder  n_genes sample tissue_ontology_term_id   
AGTACCATCCATCTGC-C0060            0   3272.0  C0060          UBERON:0013540  \
ACATCAGTCTTTAGTC-C0018            0    999.0  C0018          UBERON:0013540   
CCAATCCAGAATAGGG-C0018            0    729.0  C0018          UBERON:0013540   
CTAGCCTAGCGTTTAC-C0018            0    792.0  C0018          UBERON:0013540   
CTGAAACGTCCCTACT-C0018            0    769.0  C0018          UBERON:0013540   

                       assay_ontology_term_id cell_type_ontology_term_id   
AGTACCATCCATCTGC-C0060            EFO:0009922                 CL:0000127  \
ACATCAGTCTTTAGTC-C0018            EFO:0009899                 CL:0000127   
CCAATCCAGAATAGGG-C0018            EFO:0009899                 CL:0000127   
CTAGCCTAGCGTTTAC-C0018            EFO:0009899                 CL:0000127   
CTGAAACGTCCCTACT-C0018            EFO:0009899                 CL:0000127   

                       development_stage_ontology_term_id   
AGTACCATCCATCTGC-C0060                     HsapDv:0000165  \
ACATCAGTCTTTAGTC-C0018                     HsapDv:0000155   
CCAATCCAGAATAGGG-C0018                     HsapDv:0000155   
CTAGCCTAGCGTTTAC-C0018                     HsapDv:0000155   
CTGAAACGTCCCTACT-C0018                     HsapDv:0000155   

                       self_reported_ethnicity_ontology_term_id   
AGTACCATCCATCTGC-C0060                                  unknown  \
ACATCAGTCTTTAGTC-C0018                                  unknown   
CCAATCCAGAATAGGG-C0018                                  unknown   
CTAGCCTAGCGTTTAC-C0018                                  unknown   
CTGAAACGTCCCTACT-C0018                                  unknown   

                       disease_ontology_term_id sex_ontology_term_id   
AGTACCATCCATCTGC-C0060            MONDO:0004975         PATO:0000384  \
ACATCAGTCTTTAGTC-C0018             PATO:0000461         PATO:0000384   
CCAATCCAGAATAGGG-C0018             PATO:0000461         PATO:0000384   
CTAGCCTAGCGTTTAC-C0018             PATO:0000461         PATO:0000384   
CTGAAACGTCCCTACT-C0018             PATO:0000461         PATO:0000384   

                        is_primary_data organism_ontology_term_id   
AGTACCATCCATCTGC-C0060             True            NCBITaxon:9606  \
ACATCAGTCTTTAGTC-C0018             True            NCBITaxon:9606   
CCAATCCAGAATAGGG-C0018             True            NCBITaxon:9606   
CTAGCCTAGCGTTTAC-C0018             True            NCBITaxon:9606   
CTGAAACGTCCCTACT-C0018             True            NCBITaxon:9606   

                       suspension_type donor_id     Age_group   
AGTACCATCCATCTGC-C0060         nucleus       19      70 to 80  \
ACATCAGTCTTTAGTC-C0018         nucleus       50  less than 70   
CCAATCCAGAATAGGG-C0018         nucleus       50  less than 70   
CTAGCCTAGCGTTTAC-C0018         nucleus       50  less than 70   
CTGAAACGTCCCTACT-C0018         nucleus       50  less than 70   

                                 PMI_group RIN_group Major_celltypes   
AGTACCATCCATCTGC-C0060             unknown  4.2 to 5      Astrocytes  \
ACATCAGTCTTTAGTC-C0018  15.2 to 19.5 hours  6.1 to 7      Astrocytes   
CCAATCCAGAATAGGG-C0018  15.2 to 19.5 hours  6.1 to 7      Astrocytes   
CTAGCCTAGCGTTTAC-C0018  15.2 to 19.5 hours  6.1 to 7      Astrocytes   
CTGAAACGTCCCTACT-C0018  15.2 to 19.5 hours  6.1 to 7      Astrocytes   

                            Author_Annotation NP.Diagonis tissue_type   
AGTACCATCCATCTGC-C0060    Astrocyte-GFAP-OSMR          AD      tissue  \
ACATCAGTCTTTAGTC-C0018  Astrocyte-SLC1A2-WIF1        PART      tissue   
CCAATCCAGAATAGGG-C0018  Astrocyte-SLC1A2-SMTN        PART      tissue   
CTAGCCTAGCGTTTAC-C0018  Astrocyte-SLC1A2-WIF1        PART      tissue   
CTGAAACGTCCCTACT-C0018  Astrocyte-SLC1A2-SMTN        PART      tissue   

                        cell_type      assay            disease      organism   
AGTACCATCCATCTGC-C0060  astrocyte  10x 3' v3  Alzheimer disease  Homo sapiens  \
ACATCAGTCTTTAGTC-C0018  astrocyte  10x 3' v2             normal  Homo sapiens   
CCAATCCAGAATAGGG-C0018  astrocyte  10x 3' v2             normal  Homo sapiens   
CTAGCCTAGCGTTTAC-C0018  astrocyte  10x 3' v2             normal  Homo sapiens   
CTGAAACGTCCCTACT-C0018  astrocyte  10x 3' v2             normal  Homo sapiens   

                         sex                  tissue self_reported_ethnicity   
AGTACCATCCATCTGC-C0060  male  Brodmann (1909) area 9                 unknown  \
ACATCAGTCTTTAGTC-C0018  male  Brodmann (1909) area 9                 unknown   
CCAATCCAGAATAGGG-C0018  male  Brodmann (1909) area 9                 unknown   
CTAGCCTAGCGTTTAC-C0018  male  Brodmann (1909) area 9                 unknown   
CTGAAACGTCCCTACT-C0018  male  Brodmann (1909) area 9                 unknown   

                        development_stage observation_joinid  
AGTACCATCCATCTGC-C0060  71-year-old stage         ;H|J`@{xnY  
ACATCAGTCTTTAGTC-C0018  61-year-old stage         JQP1;<W}tH  
CCAATCCAGAATAGGG-C0018  61-year-old stage         Q3`kwi_dnf  
CTAGCCTAGCGTTTAC-C0018  61-year-old stage         QG4U-WyR#w  
CTGAAACGTCCCTACT-C0018  61-year-old stage         tcM-x)<@#d  
Size of adata.obs:
(424528, 44)

## Conclusion
- Subset for normal 
- Separate out by 3 tissues: 
	'Brodmann (1909) area 9'
	'Brodmann (1909) area 7'
	'Brodmann (1909) area 17'
- raw count is under layer `adata.raw.X.A`
- Count number of cells (normal): 85523 + 17018 + 11214 = 113,755