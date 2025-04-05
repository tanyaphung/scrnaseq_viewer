## General information
- Link: https://cellxgene.cziscience.com/collections/d0941303-7ce3-4422-9249-cf31eb98c480
- Raw counts downloaded on 2025-02-12
```
wget https://datasets.cellxgene.cziscience.com/51ed6415-8523-4bbf-9db1-c690945fa5c4.h5ad
```

## Summary
1. The downloaded h5ad file contains cells from both the Tadross study and Siletti study
- Subset for per study: 
  - 514_Tadross_Human_2025_Hypothalamus_Tadross
  - 515_Tadross_Human_2025_Hypothalamus_Siletti
2. There are 21 (or 20 because Thalamaus is misspelt and it should be Thalamus) regions of the thalamus but when subsetting for these regions, for some regions, there is only 1 cell type, which is not runnable by MAGMA gene property
- Thus, for this dataset, I'm not splitting by region
3. Cell type annotation
- Use `cell_type` as level 1: 
  - 8 cell types: 'astrocyte', 'ependymal cell', 'oligodendrocyte', 'neuron', 'endothelial cell', 'microglial cell', 'fibroblast', 'mural cell'
- Use `C2_named` as level 2: 
  - C2-1 Astrocytes-1
  - C2-4 Ependymal
  - C2-2 Astrocytes-2
  - C2-3 Astrocytes-3
  - C2-5 Oligo-Precursor-1
  - C2-7 Oligo-Mature-1
  - C2-8 Oligo-Mature-2
  - C2-6 Oligo-Precursor-2
  - C2-9 Pre-1 GABA-8
  - C2-10 Pre-1 GABA-1
  - C2-11 Pre-1 GABA-2
  - C2-12 Pre-1 GABA-3
  - C2-13 Pre-1 GABA-4
  - C2-17 Mid-1 GABA-1
  - C2-14 Pre-1 GABA-5
  - C2-25 Pre-2 GABA-1
  - C2-26 Pre-2 GABA-2
  - C2-15 Pre-1 GABA-6
  - C2-18 Mid-1 GABA-2
  - C2-27 Pre-2 GABA-3
  - C2-19 Mid-1 GABA-3
  - C2-30 Post-1 GLU-1
  - C2-16 Pre-1 GABA-7
  - C2-28 Pre-2 CHOL-1
  - C2-31 Post-1 GLU-2
  - C2-34 Mid-3 GLU-1
  - C2-29 Pre-2 GABA-4
  - C2-20 Mid-1 GABA-4
  - C2-39 Mid-2 GLU-1
  - C2-21 Mid-1 GABA-5
  - C2-40 Mid-2 GLU-2
  - C2-22 Mid-1 GABA-6
  - C2-35 Mid-3 GLU-2
  - C2-41 Mid-2 GABA-GLU-1
  - C2-36 Mid-3 GLU-3
  - C2-42 Mid-2 GABA-GLU-2
  - C2-37 Mid-3 GLU-4
  - C2-43 Mid-2 GLU-3
  - C2-44 Mid-2 GABA-GLU-3
  - C2-23 Mid-1 GABA-7
  - C2-45 Mid-2 GLU-4
  - C2-32 Post-1 GLU-3
  - C2-24 Mid-1 GABA-8
  - C2-46 Post-2 GLU-1
  - C2-47 Post-2 GLU-2
  - C2-38 Mid-3 GLU-5
  - C2-33 Post-1 GLU-4
  - C2-48 Vascular
  - C2-49 Immune Microglia
  - C2-50 Immune Macrophages
  - C2-51 Immune Myocytes
  - C2-52 Immune Tcells
  
- **NOTES:** keep in mind that there are additional cell type annotations (`C0_named`, `C1_named`, `C3_named`, and `C4_named`)

## Data exploration
### Viewing the columns of the obs.
Index(['assay_ontology_term_id', 'cell_type_ontology_term_id',
       'development_stage_ontology_term_id', 'donor_id',
       'disease_ontology_term_id', 'is_primary_data',
       'organism_ontology_term_id', 'self_reported_ethnicity_ontology_term_id',
       'sex_ontology_term_id', 'suspension_type', 'tissue_type',
       'tissue_ontology_term_id', 'nCount_RNA', 'nFeature_RNA', 'percent_mt',
       'Sample_ID', 'Dataset', 'C0_named', 'C1_named', 'C2_named', 'C3_named',
       'C4_named', 'region', 'cell_type', 'assay', 'disease', 'organism',
       'sex', 'tissue', 'self_reported_ethnicity', 'development_stage',
       'observation_joinid'],
      dtype='object')
### View the unique values.
Find unique values for column:  assay_ontology_term_id
['EFO:0009922']
Categories (1, object): ['EFO:0009922']
Find unique values for column:  cell_type_ontology_term_id
['CL:0000127', 'CL:0000065', 'CL:0000128', 'CL:0000540', 'CL:0000115', 'CL:0000129', 'CL:0000057', 'CL:0008034']
Categories (8, object): ['CL:0000057', 'CL:0000065', 'CL:0000115', 'CL:0000127', 'CL:0000128',
                         'CL:0000129', 'CL:0000540', 'CL:0008034']
Find unique values for column:  development_stage_ontology_term_id
['HsapDv:0000153', 'HsapDv:0000161', 'HsapDv:0000220', 'HsapDv:0000209', 'HsapDv:0000218', ..., 'HsapDv:0000157', 'HsapDv:0000217', 'HsapDv:0000144', 'HsapDv:0000136', 'HsapDv:0000123']
Length: 11
Categories (11, object): ['HsapDv:0000123', 'HsapDv:0000136', 'HsapDv:0000144', 'HsapDv:0000153', ...,
                          'HsapDv:0000214', 'HsapDv:0000217', 'HsapDv:0000218', 'HsapDv:0000220']
Find unique values for column:  donor_id
['TweV8', 'eXCJJ', 'tzan2', 'znZv1', 'f5sVM', ..., '3u5kk', '1C41i', 'siletti_H18.30.002', 'siletti_H19.30.001', 'siletti_H19.30.002']
Length: 11
Categories (11, object): ['1C41i', '3u5kk', 'TweV8', 'buFpQ', ..., 'siletti_H19.30.001',
                          'siletti_H19.30.002', 'tzan2', 'znZv1']
Find unique values for column:  disease_ontology_term_id
['PATO:0000461']
Categories (1, object): ['PATO:0000461']
Find unique values for column:  is_primary_data
[ True False]
Find unique values for column:  organism_ontology_term_id
['NCBITaxon:9606']
Categories (1, object): ['NCBITaxon:9606']
Find unique values for column:  self_reported_ethnicity_ontology_term_id
['unknown', 'HANCESTRO:0005']
Categories (2, object): ['unknown', 'HANCESTRO:0005']
Find unique values for column:  sex_ontology_term_id
['PATO:0000383', 'PATO:0000384']
Categories (2, object): ['PATO:0000383', 'PATO:0000384']
Find unique values for column:  suspension_type
['nucleus']
Categories (1, object): ['nucleus']
Find unique values for column:  tissue_type
['tissue']
Categories (1, object): ['tissue']
Find unique values for column:  tissue_ontology_term_id
['UBERON:0001898']
Categories (1, object): ['UBERON:0001898']
Find unique values for column:  nCount_RNA
[24137. 20008. 16589. ... 39087. 29407. 19586.]
Find unique values for column:  nFeature_RNA
[ 5853  5508  5081 ... 11617 13926 10486]
Find unique values for column:  percent_mt
[2.9829723e-01 2.8988406e-01 1.0971125e+00 ... 3.7780872e-03 2.7310925e-02
 7.9032162e-04]
Find unique values for column:  Sample_ID
['TweV8_N1', 'TweV8_US1', 'eXCJJ_N2', 'eXCJJ_US2', 'tzan2_N6.3', ..., 'siletti_10X362_5', 'siletti_10X380_5', 'siletti_10X193_4', 'siletti_10X193_3', 'siletti_10X389_2']
Length: 82
Categories (82, object): ['1C41i_A4', '1C41i_A6', '1C41i_A7', '1C41i_A8', ..., 'znZv1_S8', 'znZv1_S9',
                          'znZv1_S10', 'znZv1_S11']
Find unique values for column:  Dataset
['Tadross', 'Siletti']
Categories (2, object): ['Siletti', 'Tadross']
Find unique values for column:  C0_named
['C0-1 AstroEpendymal', 'C0-2 Oligodendrocytes', 'C0-3 Neuron', 'C0-4 Immune/Vascular']
Categories (4, object): ['C0-1 AstroEpendymal', 'C0-2 Oligodendrocytes', 'C0-3 Neuron',
                         'C0-4 Immune/Vascular']
Find unique values for column:  C1_named
['C1-1 Astrocytes', 'C1-2 Ependymal', 'C1-3 Oligo-Precursor', 'C1-4 Oligo-Mature', 'C1-5 Pre-1', ..., 'C1-9 Mid-3', 'C1-10 Mid-2', 'C1-11 Post-2', 'C1-12 Vascular', 'C1-13 Immune']
Length: 13
Categories (13, object): ['C1-1 Astrocytes', 'C1-2 Ependymal', 'C1-3 Oligo-Precursor',
                          'C1-4 Oligo-Mature', ..., 'C1-10 Mid-2', 'C1-11 Post-2', 'C1-12 Vascular',
                          'C1-13 Immune']
Find unique values for column:  C2_named
['C2-1 Astrocytes-1', 'C2-4 Ependymal', 'C2-2 Astrocytes-2', 'C2-3 Astrocytes-3', 'C2-5 Oligo-Precursor-1', ..., 'C2-48 Vascular', 'C2-49 Immune Microglia', 'C2-50 Immune Macrophages', 'C2-51 Immune Myocytes', 'C2-52 Immune Tcells']
Length: 52
Categories (52, object): ['C2-1 Astrocytes-1', 'C2-2 Astrocytes-2', 'C2-3 Astrocytes-3',
                          'C2-4 Ependymal', ..., 'C2-49 Immune Microglia',
                          'C2-50 Immune Macrophages', 'C2-51 Immune Myocytes', 'C2-52 Immune Tcells']
Find unique values for column:  C3_named
['C3-1 Astrocytes-1', 'C3-12 Ependymal Ependymocytes', 'C3-2 Astrocytes-2 VAV3', 'C3-10 Astrocytes-3 NRXN3', 'C3-3 Astrocytes-2 CPAMD8', ..., 'C3-144 Vascular Pericytes', 'C3-151 Immune Microglia SPON1', 'C3-152 Immune Microglia SLC2A3', 'C3-145 Vascular SMCs', 'C3-153 Immune Microglia STMN2']
Length: 156
Categories (156, object): ['C3-1 Astrocytes-1', 'C3-2 Astrocytes-2 VAV3', 'C3-3 Astrocytes-2 CPAMD8',
                           'C3-4 Astrocytes-2 NEFM', ..., 'C3-153 Immune Microglia STMN2',
                           'C3-154 Immune Macrophages', 'C3-155 Immune Myocytes', 'C3-156 Immune Tcells']
Find unique values for column:  C4_named
['C3-1 Astrocytes-1', 'C3-12 Ependymal Ependymocytes', 'C3-2 Astrocytes-2 VAV3', 'C3-10 Astrocytes-3 NRXN3', 'C3-3 Astrocytes-2 CPAMD8', ..., 'C3-144 Vascular Pericytes', 'C3-151 Immune Microglia SPON1', 'C3-152 Immune Microglia SLC2A3', 'C3-145 Vascular SMCs', 'C3-153 Immune Microglia STMN2']
Length: 452
Categories (452, object): ['C3-1 Astrocytes-1', 'C3-2 Astrocytes-2 VAV3', 'C3-3 Astrocytes-2 CPAMD8',
                           'C3-4 Astrocytes-2 NEFM', ..., 'C4-412 Post-2 GLU-1 SALL3 TAC3',
                           'C4-413 Post-2 GLU-1 SALL3 DMRTA2', 'C4-414 Post-2 GLU-2 TAC3',
                           'C4-415 Post-2 GLU-2 BTNL9']
Find unique values for column:  region
['NA', 'Vent', 'Fx/OT/ac', 'POA', 'Perivent', ..., 'Thalamus', 'TMN', 'Thalamaus', 'Vascular', 'ME']
Length: 21
Categories (21, object): ['ARC', 'DMH', 'Fx/OT/ac', 'LH', ..., 'Thalamus', 'VMH', 'Vascular', 'Vent']
Find unique values for column:  cell_type
['astrocyte', 'ependymal cell', 'oligodendrocyte', 'neuron', 'endothelial cell', 'microglial cell', 'fibroblast', 'mural cell']
Categories (8, object): ['fibroblast', 'ependymal cell', 'endothelial cell', 'astrocyte',
                         'oligodendrocyte', 'microglial cell', 'neuron', 'mural cell']
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
['female', 'male']
Categories (2, object): ['female', 'male']
Find unique values for column:  tissue
['hypothalamus']
Categories (1, object): ['hypothalamus']
Find unique values for column:  self_reported_ethnicity
['unknown', 'European']
Categories (2, object): ['unknown', 'European']
Find unique values for column:  development_stage
['59-year-old stage', '67-year-old stage', '94-year-old stage', '83-year-old stage', '92-year-old stage', ..., '63-year-old stage', '91-year-old stage', '50-year-old stage', '42-year-old stage', '29-year-old stage']
Length: 11
Categories (11, object): ['29-year-old stage', '42-year-old stage', '50-year-old stage',
                          '59-year-old stage', ..., '88-year-old stage', '91-year-old stage',
                          '92-year-old stage', '94-year-old stage']
Find unique values for column:  observation_joinid
['Lq%pv<Fr^H' 'fhR8Ntj$iZ' 'q|CsciwkBQ' ... 'va7u4H+!>T' 'k2~>s&>#!4'
 '%)_Pdv=(^l']
### View the adata.var:
                 feature_is_filtered feature_name feature_reference  \
ENSG00000243485                False  MIR1302-2HG    NCBITaxon:9606   
ENSG00000284662                False       OR4F16    NCBITaxon:9606   
ENSG00000237491                False    LINC01409    NCBITaxon:9606   
ENSG00000177757                False       FAM87B    NCBITaxon:9606   
ENSG00000228794                False    LINC01128    NCBITaxon:9606   

                feature_biotype feature_length    feature_type  
ENSG00000243485            gene            623          lncRNA  
ENSG00000284662            gene            939  protein_coding  
ENSG00000237491            gene           1059          lncRNA  
ENSG00000177757            gene           1947          lncRNA  
ENSG00000228794            gene           1627          lncRNA  
Size of adata.var:
(37050, 6)
### Check which layer stores the raw count:
Trying to print adata.X.A[1:25, 1:25]
[[0.         0.         0.         0.40533182 0.         0.
  0.         0.         0.         0.         0.         0.
  0.40533182 0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.47175777 0.         0.7910078  0.
  0.         0.47175777 0.         0.         0.         0.        ]
 [0.         0.5062021  0.         0.         0.         0.
  0.         0.         0.         0.5062021  0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.86497575 0.         1.1192063  0.
  0.         0.         0.         0.         0.5232329  0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.5401897  0.         0.         0.
  0.5401897  0.5401897  0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.9395091  0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.9441867  0.         0.         0.
  0.         0.         0.         0.         0.9441867  0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.588784   0.        ]
 [0.         0.5899997  0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.5899997  0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         1.2708508  0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  1.0497154  0.         0.65670043 0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.6688174  0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         1.0660393  0.
  0.         0.         0.         0.         0.6688174  0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         1.3588834  0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  1.0766383  0.         0.         0.         0.6767122  0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.711014   0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.730035   0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.7463383  0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         1.1904594  0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.79615    0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.840499   0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.840499   0.        ]
 [0.         0.88278997 0.         0.88278997 0.         0.
  0.88278997 0.         0.         0.         0.         0.
  0.88278997 0.         0.         0.         0.88278997 0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.89604956 0.         0.
  0.         0.         0.89604956 0.         0.         0.
  1.6770437  0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.9154515  0.         0.         0.
  0.9154515  0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]]
Trying to print adata.X[1:25, 1:25]
  (0, 3)	0.40533182
  (0, 12)	0.40533182
  (1, 14)	0.47175777
  (1, 16)	0.7910078
  (1, 19)	0.47175777
  (2, 1)	0.5062021
  (2, 9)	0.5062021
  (3, 14)	0.86497575
  (3, 16)	1.1192063
  (3, 22)	0.5232329
  (4, 8)	0.5401897
  (4, 12)	0.5401897
  (4, 13)	0.5401897
  (5, 12)	0.9395091
  (6, 8)	0.9441867
  (6, 16)	0.9441867
  (7, 22)	0.588784
  (8, 1)	0.5899997
  (8, 16)	0.5899997
  (9, 19)	1.2708508
  (10, 12)	1.0497154
  (10, 14)	0.65670043
  (11, 1)	0.6688174
  (11, 16)	1.0660393
  (11, 22)	0.6688174
  (12, 16)	1.3588834
  (13, 12)	1.0766383
  (13, 16)	0.6767122
  (14, 22)	0.711014
  (15, 12)	0.730035
  (16, 16)	0.7463383
  (17, 1)	1.1904594
  (19, 8)	0.79615
  (20, 3)	0.840499
  (20, 22)	0.840499
  (21, 1)	0.88278997
  (21, 3)	0.88278997
  (21, 6)	0.88278997
  (21, 12)	0.88278997
  (21, 16)	0.88278997
  (22, 3)	0.89604956
  (22, 8)	0.89604956
  (22, 12)	1.6770437
  (23, 8)	0.9154515
  (23, 12)	0.9154515
Trying to print adata.raw.X.A[1:25, 1:25]
[[0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 2. 0. 0. 1. 0. 0. 0. 0.]
 [0. 1. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 3. 0. 0. 0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 3. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 3. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
 [0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 1. 0. 1. 0. 0. 1. 0. 0. 0. 0. 0. 1. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 1. 0. 0. 0. 0. 1. 0. 0. 0. 3. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]]
### View the adata.obs:
                            assay_ontology_term_id cell_type_ontology_term_id  \
TweV8_N1_TGACAGTTCTCCGAGG-1            EFO:0009922                 CL:0000127   
TweV8_N1_CAGATTGCACTTTATC-1            EFO:0009922                 CL:0000127   
TweV8_N1_TCAGTTTGTTGGCTAT-1            EFO:0009922                 CL:0000065   
TweV8_N1_CTTACCGTCCTGATAG-1            EFO:0009922                 CL:0000127   
TweV8_N1_AGACAAAAGAGCTGCA-1            EFO:0009922                 CL:0000127   

                            development_stage_ontology_term_id donor_id  \
TweV8_N1_TGACAGTTCTCCGAGG-1                     HsapDv:0000153    TweV8   
TweV8_N1_CAGATTGCACTTTATC-1                     HsapDv:0000153    TweV8   
TweV8_N1_TCAGTTTGTTGGCTAT-1                     HsapDv:0000153    TweV8   
TweV8_N1_CTTACCGTCCTGATAG-1                     HsapDv:0000153    TweV8   
TweV8_N1_AGACAAAAGAGCTGCA-1                     HsapDv:0000153    TweV8   

                            disease_ontology_term_id  is_primary_data  \
TweV8_N1_TGACAGTTCTCCGAGG-1             PATO:0000461             True   
TweV8_N1_CAGATTGCACTTTATC-1             PATO:0000461             True   
TweV8_N1_TCAGTTTGTTGGCTAT-1             PATO:0000461             True   
TweV8_N1_CTTACCGTCCTGATAG-1             PATO:0000461             True   
TweV8_N1_AGACAAAAGAGCTGCA-1             PATO:0000461             True   

                            organism_ontology_term_id  \
TweV8_N1_TGACAGTTCTCCGAGG-1            NCBITaxon:9606   
TweV8_N1_CAGATTGCACTTTATC-1            NCBITaxon:9606   
TweV8_N1_TCAGTTTGTTGGCTAT-1            NCBITaxon:9606   
TweV8_N1_CTTACCGTCCTGATAG-1            NCBITaxon:9606   
TweV8_N1_AGACAAAAGAGCTGCA-1            NCBITaxon:9606   

                            self_reported_ethnicity_ontology_term_id  \
TweV8_N1_TGACAGTTCTCCGAGG-1                                  unknown   
TweV8_N1_CAGATTGCACTTTATC-1                                  unknown   
TweV8_N1_TCAGTTTGTTGGCTAT-1                                  unknown   
TweV8_N1_CTTACCGTCCTGATAG-1                                  unknown   
TweV8_N1_AGACAAAAGAGCTGCA-1                                  unknown   

                            sex_ontology_term_id suspension_type tissue_type  \
TweV8_N1_TGACAGTTCTCCGAGG-1         PATO:0000383         nucleus      tissue   
TweV8_N1_CAGATTGCACTTTATC-1         PATO:0000383         nucleus      tissue   
TweV8_N1_TCAGTTTGTTGGCTAT-1         PATO:0000383         nucleus      tissue   
TweV8_N1_CTTACCGTCCTGATAG-1         PATO:0000383         nucleus      tissue   
TweV8_N1_AGACAAAAGAGCTGCA-1         PATO:0000383         nucleus      tissue   

                            tissue_ontology_term_id  nCount_RNA  nFeature_RNA  \
TweV8_N1_TGACAGTTCTCCGAGG-1          UBERON:0001898     24137.0          5853   
TweV8_N1_CAGATTGCACTTTATC-1          UBERON:0001898     20008.0          5508   
TweV8_N1_TCAGTTTGTTGGCTAT-1          UBERON:0001898     16589.0          5081   
TweV8_N1_CTTACCGTCCTGATAG-1          UBERON:0001898     15175.0          4613   
TweV8_N1_AGACAAAAGAGCTGCA-1          UBERON:0001898     14546.0          4802   

                             percent_mt Sample_ID  Dataset  \
TweV8_N1_TGACAGTTCTCCGAGG-1    0.298297  TweV8_N1  Tadross   
TweV8_N1_CAGATTGCACTTTATC-1    0.289884  TweV8_N1  Tadross   
TweV8_N1_TCAGTTTGTTGGCTAT-1    1.097113  TweV8_N1  Tadross   
TweV8_N1_CTTACCGTCCTGATAG-1    0.019769  TweV8_N1  Tadross   
TweV8_N1_AGACAAAAGAGCTGCA-1    0.219992  TweV8_N1  Tadross   

                                        C0_named         C1_named  \
TweV8_N1_TGACAGTTCTCCGAGG-1  C0-1 AstroEpendymal  C1-1 Astrocytes   
TweV8_N1_CAGATTGCACTTTATC-1  C0-1 AstroEpendymal  C1-1 Astrocytes   
TweV8_N1_TCAGTTTGTTGGCTAT-1  C0-1 AstroEpendymal   C1-2 Ependymal   
TweV8_N1_CTTACCGTCCTGATAG-1  C0-1 AstroEpendymal  C1-1 Astrocytes   
TweV8_N1_AGACAAAAGAGCTGCA-1  C0-1 AstroEpendymal  C1-1 Astrocytes   

                                      C2_named                       C3_named  \
TweV8_N1_TGACAGTTCTCCGAGG-1  C2-1 Astrocytes-1              C3-1 Astrocytes-1   
TweV8_N1_CAGATTGCACTTTATC-1  C2-1 Astrocytes-1              C3-1 Astrocytes-1   
TweV8_N1_TCAGTTTGTTGGCTAT-1     C2-4 Ependymal  C3-12 Ependymal Ependymocytes   
TweV8_N1_CTTACCGTCCTGATAG-1  C2-1 Astrocytes-1              C3-1 Astrocytes-1   
TweV8_N1_AGACAAAAGAGCTGCA-1  C2-2 Astrocytes-2         C3-2 Astrocytes-2 VAV3   

                                                  C4_named region  \
TweV8_N1_TGACAGTTCTCCGAGG-1              C3-1 Astrocytes-1     NA   
TweV8_N1_CAGATTGCACTTTATC-1              C3-1 Astrocytes-1     NA   
TweV8_N1_TCAGTTTGTTGGCTAT-1  C3-12 Ependymal Ependymocytes   Vent   
TweV8_N1_CTTACCGTCCTGATAG-1              C3-1 Astrocytes-1     NA   
TweV8_N1_AGACAAAAGAGCTGCA-1         C3-2 Astrocytes-2 VAV3     NA   

                                  cell_type      assay disease      organism  \
TweV8_N1_TGACAGTTCTCCGAGG-1       astrocyte  10x 3' v3  normal  Homo sapiens   
TweV8_N1_CAGATTGCACTTTATC-1       astrocyte  10x 3' v3  normal  Homo sapiens   
TweV8_N1_TCAGTTTGTTGGCTAT-1  ependymal cell  10x 3' v3  normal  Homo sapiens   
TweV8_N1_CTTACCGTCCTGATAG-1       astrocyte  10x 3' v3  normal  Homo sapiens   
TweV8_N1_AGACAAAAGAGCTGCA-1       astrocyte  10x 3' v3  normal  Homo sapiens   

                                sex        tissue self_reported_ethnicity  \
TweV8_N1_TGACAGTTCTCCGAGG-1  female  hypothalamus                 unknown   
TweV8_N1_CAGATTGCACTTTATC-1  female  hypothalamus                 unknown   
TweV8_N1_TCAGTTTGTTGGCTAT-1  female  hypothalamus                 unknown   
TweV8_N1_CTTACCGTCCTGATAG-1  female  hypothalamus                 unknown   
TweV8_N1_AGACAAAAGAGCTGCA-1  female  hypothalamus                 unknown   

                             development_stage observation_joinid  
TweV8_N1_TGACAGTTCTCCGAGG-1  59-year-old stage         Lq%pv<Fr^H  
TweV8_N1_CAGATTGCACTTTATC-1  59-year-old stage         fhR8Ntj$iZ  
TweV8_N1_TCAGTTTGTTGGCTAT-1  59-year-old stage         q|CsciwkBQ  
TweV8_N1_CTTACCGTCCTGATAG-1  59-year-old stage         XJ4Rzc1$-O  
TweV8_N1_AGACAAAAGAGCTGCA-1  59-year-old stage         2&d2s@3ydN  
Size of adata.obs:
(433369, 32)