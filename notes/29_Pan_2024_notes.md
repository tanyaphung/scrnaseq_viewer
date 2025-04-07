## General information
- Link: https://cellxgene.cziscience.com/collections/7c4552fd-8a6d-4da3-9854-2dfa8baca8bf
- Raw counts downloaded on 2025-02-13

```
wget https://datasets.cellxgene.cziscience.com/982952fb-dee0-4498-a072-3dc498b06aae.h5ad
```

### Viewing the columns of the obs.
Index(['APOE', 'Brain.Region', 'SORT', 'Braak.stage', 'Amyloid',
       'Brain.weight', 'PMI.hr.', 'RIN', 'Gray.vs.White', 'n_genes_by_counts',
       'total_counts', 'pct_counts_ribo', 'pct_counts_hb', 'pct_counts_mt',
       'n_counts', 'n_genes', '20_leiden_1.0', 'Sample', 'Molecule',
       'single_or_paired_end', 'Instrument', 'organism_ontology_term_id',
       'donor_id', 'development_stage_ontology_term_id',
       'sex_ontology_term_id', 'self_reported_ethnicity_ontology_term_id',
       'disease_ontology_term_id', 'tissue_type', 'tissue_ontology_term_id',
       'cell_type_ontology_term_id', 'assay_ontology_term_id',
       'suspension_type', 'is_primary_data', 'cell_type', 'assay', 'disease',
       'organism', 'sex', 'tissue', 'self_reported_ethnicity',
       'development_stage', 'observation_joinid'],
      dtype='object')
### View the unique values.
Find unique values for column:  APOE
['E3/E4', 'E3/E3', 'E2/E3']
Categories (3, object): ['E2/E3', 'E3/E3', 'E3/E4']
Find unique values for column:  Brain.Region
['Medial temporal', 'Occipital', 'Parietal', 'Frontal']
Categories (4, object): ['Medial temporal', 'Occipital', 'Parietal', 'Frontal']
Find unique values for column:  SORT
['DAPI', 'DAPI/NeuNneg', 'Pax6/NeuNneg']
Categories (3, object): ['DAPI', 'DAPI/NeuNneg', 'Pax6/NeuNneg']
Find unique values for column:  Braak.stage
['IV', '0', 'III', 'I']
Categories (4, object): ['0', 'I', 'III', 'IV']
Find unique values for column:  Amyloid
['C1', 'C0', 'Diffuse plaques']
Categories (3, object): ['C1', 'C0', 'Diffuse plaques']
Find unique values for column:  Brain.weight
[1355.   nan 1280. 1375.]
Find unique values for column:  PMI.hr.
[24.   15.2  15.   19.5  11.   13.    3.75 10.  ]
Find unique values for column:  RIN
['4.3', '4.1', '6.2', '4.6', '4.5', ..., '5.2', '5.4', '5.3', '5.7', '6.1']
Length: 12
Categories (12, object): ['5', '5.2', '5.3', '5.4', ..., '6.2', '4.6', '4.5', '6.1']
Find unique values for column:  Gray.vs.White
['wm', 'gray', 'gray/wm']
Categories (3, object): ['gray', 'gray/wm', 'wm']
Find unique values for column:  n_genes_by_counts
[ 991  763  581 ... 6209 8427 3908]
Find unique values for column:  total_counts
[ 1764  1251   848 ... 12557 16714 15810]
Find unique values for column:  pct_counts_ribo
[0.45351475 0.23980816 0.23584905 ... 0.36875945 0.34537342 0.50257266]
Find unique values for column:  pct_counts_hb
[0.         0.07993605 0.03016591 ... 0.03940887 0.00598301 0.00632511]
Find unique values for column:  pct_counts_mt
[1.3038548  5.8353314  1.0613208  ... 0.11003791 1.0193517  0.4247936 ]
Find unique values for column:  n_counts
[ 1764.  1251.   848. ... 12557. 16714. 15810.]
Find unique values for column:  n_genes
[ 991  763  581 ... 6209 8427 3908]
Find unique values for column:  20_leiden_1.0
['0', '7', '5', '10', '9', ..., '23', '24', '21', '22', '25']
Length: 26
Categories (26, object): ['0', '1', '2', '3', ..., '22', '23', '24', '25']
Find unique values for column:  Sample
['ALSP MTL wm', 'ALSP OCC wm', 'ALSP PAR wm', 'ALSP PF gm', 'ALSP PF wm', ..., 'AD 1', 'Control 1', 'Control 2', 'Control 3', 'AD 2']
Length: 12
Categories (12, object): ['ALSP PF gm' < 'ALSP PF wm' < 'ALSP MTL wm' < 'ALSP OCC wm' ... 'Control 4' <
                          'Control 5' < 'AD 1' < 'AD 2']
Find unique values for column:  Molecule
['mRNA']
Categories (1, object): ['mRNA']
Find unique values for column:  single_or_paired_end
['paired end']
Categories (1, object): ['paired end']
Find unique values for column:  Instrument
['Novaseq6000 S4 2x150']
Categories (1, object): ['Novaseq6000 S4 2x150']
Find unique values for column:  organism_ontology_term_id
['NCBITaxon:9606']
Categories (1, object): ['NCBITaxon:9606']
Find unique values for column:  donor_id
['ALSP', '108', '100', '46', '94', '112', '116', '30']
Categories (8, object): ['46', '94', '100', '108', '112', '116', '30', 'ALSP']
Find unique values for column:  development_stage_ontology_term_id
['HsapDv:0000166', 'HsapDv:0000170', 'HsapDv:0000147', 'HsapDv:0000209', 'HsapDv:0000208', 'HsapDv:0000164', 'HsapDv:0000162', 'HsapDv:0000207']
Categories (8, object): ['HsapDv:0000147', 'HsapDv:0000162', 'HsapDv:0000164', 'HsapDv:0000166',
                         'HsapDv:0000170', 'HsapDv:0000208', 'HsapDv:0000209', 'HsapDv:0000207']
Find unique values for column:  sex_ontology_term_id
['PATO:0000384', 'PATO:0000383']
Categories (2, object): ['PATO:0000383', 'PATO:0000384']
Find unique values for column:  self_reported_ethnicity_ontology_term_id
['HANCESTRO:0005', 'HANCESTRO:0016', 'unknown']
Categories (3, object): ['HANCESTRO:0005', 'HANCESTRO:0016', 'unknown']
Find unique values for column:  disease_ontology_term_id
['MONDO:0800027', 'PATO:0000461', 'MONDO:0004975']
Categories (3, object): ['MONDO:0004975', 'MONDO:0800027', 'PATO:0000461']
Find unique values for column:  tissue_type
['tissue']
Categories (1, object): ['tissue']
Find unique values for column:  tissue_ontology_term_id
['UBERON:0016534', 'UBERON:0016535', 'UBERON:0016531', 'UBERON:0000451', 'UBERON:0016528']
Categories (5, object): ['UBERON:0000451', 'UBERON:0016535', 'UBERON:0016534', 'UBERON:0016531',
                         'UBERON:0016528']
Find unique values for column:  cell_type_ontology_term_id
['CL:0000128', 'CL:0000127', 'CL:0000129', 'CL:0000115', 'CL:4023051', 'CL:0002453', 'CL:0000498', 'CL:0000679']
Categories (8, object): ['CL:0000115', 'CL:0000127', 'CL:0000128', 'CL:0000129', 'CL:0000498',
                         'CL:0002453', 'CL:0000679', 'CL:4023051']
Find unique values for column:  assay_ontology_term_id
['EFO:0009922']
Categories (1, object): ['EFO:0009922']
Find unique values for column:  suspension_type
['nucleus']
Categories (1, object): ['nucleus']
Find unique values for column:  is_primary_data
[ True]
Find unique values for column:  cell_type
['oligodendrocyte', 'astrocyte', 'microglial cell', 'endothelial cell', 'vascular leptomeningeal cell', 'oligodendrocyte precursor cell', 'inhibitory interneuron', 'glutamatergic neuron']
Categories (8, object): ['endothelial cell', 'astrocyte', 'oligodendrocyte', 'microglial cell',
                         'inhibitory interneuron', 'oligodendrocyte precursor cell',
                         'glutamatergic neuron', 'vascular leptomeningeal cell']
Find unique values for column:  assay
['10x 3' v3']
Categories (1, object): ['10x 3' v3']
Find unique values for column:  disease
['leukoencephalopathy, diffuse hereditary, with..., 'normal', 'Alzheimer disease']
Categories (3, object): ['Alzheimer disease', 'leukoencephalopathy, diffuse hereditary, with...,
                         'normal']
Find unique values for column:  organism
['Homo sapiens']
Categories (1, object): ['Homo sapiens']
Find unique values for column:  sex
['male', 'female']
Categories (2, object): ['female', 'male']
Find unique values for column:  tissue
['white matter of temporal lobe', 'white matter of occipital lobe', 'white matter of parietal lobe', 'prefrontal cortex', 'white matter of frontal lobe']
Categories (5, object): ['prefrontal cortex', 'white matter of occipital lobe',
                         'white matter of temporal lobe', 'white matter of parietal lobe',
                         'white matter of frontal lobe']
Find unique values for column:  self_reported_ethnicity
['European', 'African American or Afro-Caribbean', 'unknown']
Categories (3, object): ['European', 'African American or Afro-Caribbean', 'unknown']
Find unique values for column:  development_stage
['72-year-old stage', '76-year-old stage', '53-year-old stage', '83-year-old stage', '82-year-old stage', '70-year-old stage', '68-year-old stage', '81-year-old stage']
Categories (8, object): ['53-year-old stage', '68-year-old stage', '70-year-old stage',
                         '72-year-old stage', '76-year-old stage', '82-year-old stage',
                         '83-year-old stage', '81-year-old stage']
Find unique values for column:  observation_joinid
['0b*^?3z<35' 'rRUZAD#H$=' ';E{vC#|}Wl' ... 'VK#cjj5P$l' 'R8HyT$M%XD'
 '#y*WcjW%}E']
### View the adata.var:
                    mt   ribo     hb  highly_variable     means  dispersions  \
new_index                                                                      
ENSG00000186827  False  False  False            False  0.002057     2.044141   
ENSG00000186891  False  False  False            False  0.003529     1.850113   
ENSG00000160072  False  False  False             True  0.171040     1.959340   
ENSG00000260179  False  False  False            False  0.001390     1.191500   
ENSG00000234396  False  False  False            False  0.007698     1.369916   

                 dispersions_norm          gene_name  feature_is_filtered  \
new_index                                                                   
ENSG00000186827          0.651523            TNFRSF4                False   
ENSG00000186891          0.415619           TNFRSF18                False   
ENSG00000160072          0.548420             ATAD3B                False   
ENSG00000260179         -0.385141  ENSG00000260179.1                False   
ENSG00000234396         -0.168218  ENSG00000234396.3                False   

                       gene_version       feature_name feature_reference  \
new_index                                                                  
ENSG00000186827  ENSG00000186827.11            TNFRSF4    NCBITaxon:9606   
ENSG00000186891  ENSG00000186891.14           TNFRSF18    NCBITaxon:9606   
ENSG00000160072  ENSG00000160072.20             ATAD3B    NCBITaxon:9606   
ENSG00000260179   ENSG00000260179.1  ENSG00000260179.1    NCBITaxon:9606   
ENSG00000234396   ENSG00000234396.3  ENSG00000234396.3    NCBITaxon:9606   

                feature_biotype feature_length    feature_type  
new_index                                                       
ENSG00000186827            gene           1039  protein_coding  
ENSG00000186891            gene            789  protein_coding  
ENSG00000160072            gene           3300  protein_coding  
ENSG00000260179            gene           1558          lncRNA  
ENSG00000234396            gene            326          lncRNA  
Size of adata.var:
(61427, 15)
### Check which layer stores the raw count:
Trying to print adata.X.A[1:25, 1:25]
[[0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        2.5488553 0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        2.1449342 0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.9432759 0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        1.7930106 0.        0.        0.        0.        0.
  0.        0.        0.        0.        1.7930106 0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        1.3025895 0.       ]
 [0.        0.6654921 0.        0.        0.        0.        0.
  0.        0.        1.0615662 0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.6654921 0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        3.753418  0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        1.0019338 0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        2.1831274 0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        1.8116157 0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        1.7732519 0.        0.        0.        0.        0.
  0.        0.        2.7517433 0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        1.5772703 0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]
 [0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.        0.        0.        0.        0.
  0.        0.        0.       ]]
Trying to print adata.X[1:25, 1:25]
  (1, 1)	2.5488553
  (2, 9)	2.1449342
  (4, 11)	0.9432759
  (6, 1)	1.7930106
  (6, 11)	1.7930106
  (7, 22)	1.3025895
  (8, 1)	0.6654921
  (8, 9)	1.0615662
  (8, 22)	0.6654921
  (11, 11)	3.753418
  (12, 9)	1.0019338
  (13, 9)	2.1831274
  (14, 9)	1.8116157
  (16, 1)	1.7732519
  (16, 9)	2.7517433
  (22, 11)	1.5772703
Trying to print adata.raw.X.A[1:25, 1:25]
[[0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 5. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 1. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 1. 0. 0. 0. 0. 0. 0. 0. 3. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]]
### View the adata.obs:
             APOE     Brain.Region  SORT Braak.stage Amyloid  Brain.weight  \
6-MTL_WM    E3/E4  Medial temporal  DAPI          IV      C1        1355.0   
44-MTL_WM   E3/E4  Medial temporal  DAPI          IV      C1        1355.0   
77-MTL_WM   E3/E4  Medial temporal  DAPI          IV      C1        1355.0   
82-MTL_WM   E3/E4  Medial temporal  DAPI          IV      C1        1355.0   
106-MTL_WM  E3/E4  Medial temporal  DAPI          IV      C1        1355.0   

            PMI.hr.  RIN Gray.vs.White  n_genes_by_counts  total_counts  \
6-MTL_WM       24.0  4.3            wm                991          1764   
44-MTL_WM      24.0  4.3            wm                763          1251   
77-MTL_WM      24.0  4.3            wm                581           848   
82-MTL_WM      24.0  4.3            wm               2711          6630   
106-MTL_WM     24.0  4.3            wm                779          1293   

            pct_counts_ribo  pct_counts_hb  pct_counts_mt  n_counts  n_genes  \
6-MTL_WM           0.453515       0.000000       1.303855    1764.0      991   
44-MTL_WM          0.239808       0.079936       5.835331    1251.0      763   
77-MTL_WM          0.235849       0.000000       1.061321     848.0      581   
82-MTL_WM          0.316742       0.030166       1.990950    6630.0     2711   
106-MTL_WM         0.464037       0.000000       1.546790    1293.0      779   

           20_leiden_1.0       Sample Molecule single_or_paired_end  \
6-MTL_WM               0  ALSP MTL wm     mRNA           paired end   
44-MTL_WM              7  ALSP MTL wm     mRNA           paired end   
77-MTL_WM              7  ALSP MTL wm     mRNA           paired end   
82-MTL_WM              5  ALSP MTL wm     mRNA           paired end   
106-MTL_WM             7  ALSP MTL wm     mRNA           paired end   

                      Instrument organism_ontology_term_id donor_id  \
6-MTL_WM    Novaseq6000 S4 2x150            NCBITaxon:9606     ALSP   
44-MTL_WM   Novaseq6000 S4 2x150            NCBITaxon:9606     ALSP   
77-MTL_WM   Novaseq6000 S4 2x150            NCBITaxon:9606     ALSP   
82-MTL_WM   Novaseq6000 S4 2x150            NCBITaxon:9606     ALSP   
106-MTL_WM  Novaseq6000 S4 2x150            NCBITaxon:9606     ALSP   

           development_stage_ontology_term_id sex_ontology_term_id  \
6-MTL_WM                       HsapDv:0000166         PATO:0000384   
44-MTL_WM                      HsapDv:0000166         PATO:0000384   
77-MTL_WM                      HsapDv:0000166         PATO:0000384   
82-MTL_WM                      HsapDv:0000166         PATO:0000384   
106-MTL_WM                     HsapDv:0000166         PATO:0000384   

           self_reported_ethnicity_ontology_term_id disease_ontology_term_id  \
6-MTL_WM                             HANCESTRO:0005            MONDO:0800027   
44-MTL_WM                            HANCESTRO:0005            MONDO:0800027   
77-MTL_WM                            HANCESTRO:0005            MONDO:0800027   
82-MTL_WM                            HANCESTRO:0005            MONDO:0800027   
106-MTL_WM                           HANCESTRO:0005            MONDO:0800027   

           tissue_type tissue_ontology_term_id cell_type_ontology_term_id  \
6-MTL_WM        tissue          UBERON:0016534                 CL:0000128   
44-MTL_WM       tissue          UBERON:0016534                 CL:0000128   
77-MTL_WM       tissue          UBERON:0016534                 CL:0000128   
82-MTL_WM       tissue          UBERON:0016534                 CL:0000127   
106-MTL_WM      tissue          UBERON:0016534                 CL:0000128   

           assay_ontology_term_id suspension_type  is_primary_data  \
6-MTL_WM              EFO:0009922         nucleus             True   
44-MTL_WM             EFO:0009922         nucleus             True   
77-MTL_WM             EFO:0009922         nucleus             True   
82-MTL_WM             EFO:0009922         nucleus             True   
106-MTL_WM            EFO:0009922         nucleus             True   

                  cell_type      assay  \
6-MTL_WM    oligodendrocyte  10x 3' v3   
44-MTL_WM   oligodendrocyte  10x 3' v3   
77-MTL_WM   oligodendrocyte  10x 3' v3   
82-MTL_WM         astrocyte  10x 3' v3   
106-MTL_WM  oligodendrocyte  10x 3' v3   

                                                      disease      organism  \
6-MTL_WM    leukoencephalopathy, diffuse hereditary, with ...  Homo sapiens   
44-MTL_WM   leukoencephalopathy, diffuse hereditary, with ...  Homo sapiens   
77-MTL_WM   leukoencephalopathy, diffuse hereditary, with ...  Homo sapiens   
82-MTL_WM   leukoencephalopathy, diffuse hereditary, with ...  Homo sapiens   
106-MTL_WM  leukoencephalopathy, diffuse hereditary, with ...  Homo sapiens   

             sex                         tissue self_reported_ethnicity  \
6-MTL_WM    male  white matter of temporal lobe                European   
44-MTL_WM   male  white matter of temporal lobe                European   
77-MTL_WM   male  white matter of temporal lobe                European   
82-MTL_WM   male  white matter of temporal lobe                European   
106-MTL_WM  male  white matter of temporal lobe                European   

            development_stage observation_joinid  
6-MTL_WM    72-year-old stage         0b*^?3z<35  
44-MTL_WM   72-year-old stage         rRUZAD#H$=  
77-MTL_WM   72-year-old stage         ;E{vC#|}Wl  
82-MTL_WM   72-year-old stage         vYCKi%e9%T  
106-MTL_WM  72-year-old stage         r<T9CWQrcN  
Size of adata.obs:
(61747, 42)
### Number of female donors: 
3
### Number of male donors: 
5

```
adata_normal_subset = adata[cell_id_normal]
adata_disease_subset = adata[cell_id_disease]
print("Subset for normal:")
print(adata_normal_subset.obs.shape)
for i in adata_normal_subset.obs["tissue"].unique():
    print(i)


print("Subset for disease:")
print(adata_disease_subset.obs.shape)
```

```
Subset for normal:
(31971, 42)
prefrontal cortex
Subset for disease:
(29776, 42)
```
- Therefore only prefrontal cortex for normal samples