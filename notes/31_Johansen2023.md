# Documentation
## General information
- Link: https://cellxgene.cziscience.com/collections/35928d1c-36fc-4f93-9a8d-0b921ab41745
- Raw counts downloaded on 2025-03-12
    - Note that there is a file for each cell type so I will need to combine these
    - In total there are 379330 (~400k) cells. It would be computationally intensive to combine all of the h5ad into one h5ad files to have ~400k cells. Therefore, what I will be doing is to preprocess each one. Then, I will combine the files after computing mean.

```
#L23IT
wget https://datasets.cellxgene.cziscience.com/ed71fbbd-1581-4999-972d-b696d8473825.h5ad
mv ed71fbbd-1581-4999-972d-b696d8473825.h5ad L23IT.h5ad

#Oligo
wget https://datasets.cellxgene.cziscience.com/c549802b-9e3b-4c75-bdfc-d790cc322129.h5ad
mv c549802b-9e3b-4c75-bdfc-d790cc322129.h5ad Oligo.h5ad

#L4IT
wget https://datasets.cellxgene.cziscience.com/df21d007-9b08-411e-9dba-978ab4d4cb96.h5ad
mv df21d007-9b08-411e-9dba-978ab4d4cb96.h5ad L4IT.h5ad

#L5IT
wget https://datasets.cellxgene.cziscience.com/e1b881b3-8c5d-407c-abcc-add006166b4c.h5ad
mv e1b881b3-8c5d-407c-abcc-add006166b4c.h5ad L5IT.h5ad

#Sst
wget https://datasets.cellxgene.cziscience.com/a6cca164-4915-46e7-b9de-cdc348d25f3d.h5ad
mv a6cca164-4915-46e7-b9de-cdc348d25f3d.h5ad Sst.h5ad

#Vip
wget https://datasets.cellxgene.cziscience.com/3a16df1a-df06-4b7a-9f4f-66eecd71118b.h5ad
mv 3a16df1a-df06-4b7a-9f4f-66eecd71118b.h5ad Vip.h5ad

#Pvalb
wget https://datasets.cellxgene.cziscience.com/38bb351a-0e6f-4c99-9be9-80eec49c7491.h5ad
mv 38bb351a-0e6f-4c99-9be9-80eec49c7491.h5ad Pvalb.h5ad

#L6IT
wget https://datasets.cellxgene.cziscience.com/8734931c-fd2c-4af3-9d89-32a8cb4ce7a4.h5ad
mv 8734931c-fd2c-4af3-9d89-32a8cb4ce7a4.h5ad L6IT.h5ad

#Astro
wget https://datasets.cellxgene.cziscience.com/1ca450e2-f066-491a-b549-dab81c60bd10.h5ad
mv 1ca450e2-f066-491a-b549-dab81c60bd10.h5ad Astro.h5ad

#OPC
wget https://datasets.cellxgene.cziscience.com/310fffda-1c59-4880-83c3-a10861834c54.h5ad
mv 310fffda-1c59-4880-83c3-a10861834c54.h5ad OPC.h5ad

#L6b
wget https://datasets.cellxgene.cziscience.com/41545de3-4c25-498c-bd45-b0f064289269.h5ad
mv 41545de3-4c25-498c-bd45-b0f064289269.h5ad L6b.h5ad

#Lamp5
wget https://datasets.cellxgene.cziscience.com/449a9243-4a55-4334-b971-ed682755b1cc.h5ad 
mv 449a9243-4a55-4334-b971-ed682755b1cc.h5ad Lamp5.h5ad

#L6CT
wget https://datasets.cellxgene.cziscience.com/f77c7104-d924-4be0-96ad-3b410773ae24.h5ad 
mv f77c7104-d924-4be0-96ad-3b410773ae24.h5ad 

#L56NP
wget https://datasets.cellxgene.cziscience.com/ede11375-4da3-40ff-b346-54a2329a6de4.h5ad
mv ede11375-4da3-40ff-b346-54a2329a6de4.h5ad L56NP.h5ad

#L6ITCar3
wget https://datasets.cellxgene.cziscience.com/177a2189-1bab-42ec-8ff3-2d0d8b28aa44.h5ad
mv 177a2189-1bab-42ec-8ff3-2d0d8b28aa44.h5ad L6ITCar3.h5ad

#Lamp5Lhx6
wget https://datasets.cellxgene.cziscience.com/952877cd-1125-41bd-b2e5-610f187b9097.h5ad
mv 952877cd-1125-41bd-b2e5-610f187b9097.h5ad Lamp5Lhx6.h5ad

#Sncg
wget https://datasets.cellxgene.cziscience.com/c23ea1a6-b2c8-45b7-b34c-b85ea6e85a51.h5ad
mv c23ea1a6-b2c8-45b7-b34c-b85ea6e85a51.h5ad Sncg.h5ad

#MicroPVM
wget https://datasets.cellxgene.cziscience.com/fd7713fc-aac4-4dbd-853b-11b367714efa.h5ad
mv fd7713fc-aac4-4dbd-853b-11b367714efa.h5ad MicroPVM.h5ad

#Chandelier
wget https://datasets.cellxgene.cziscience.com/d6845f51-ddd0-4b90-8c43-4f1063853388.h5ad
mv d6845f51-ddd0-4b90-8c43-4f1063853388.h5ad Chandelier.h5ad

#Pax6
wget https://datasets.cellxgene.cziscience.com/bee7de3b-2e4e-488d-b091-c80be5e2d91e.h5ad
mv bee7de3b-2e4e-488d-b091-c80be5e2d91e.h5ad Pax6.h5ad

#VLMC
wget https://datasets.cellxgene.cziscience.com/b72ee73c-3824-406f-956a-92f861f4cc32.h5ad
mv b72ee73c-3824-406f-956a-92f861f4cc32.h5ad VLMC.h5ad

#Endo
wget https://datasets.cellxgene.cziscience.com/d93765a2-c22a-45fc-936f-c85d2a6865ed.h5ad
mv d93765a2-c22a-45fc-936f-c85d2a6865ed.h5ad Endo.h5ad

#L5ET
wget https://datasets.cellxgene.cziscience.com/a6156273-d49c-45d7-8b68-e2bd17592c6c.h5ad
mv a6156273-d49c-45d7-8b68-e2bd17592c6c.h5ad L5ET.h5ad

#SstChodl
wget https://datasets.cellxgene.cziscience.com/1bed4160-d84f-4c07-bad2-0a8abf7ffa53.h5ad
mv 1bed4160-d84f-4c07-bad2-0a8abf7ffa53.h5ad SstChodl.h5ad
```

### Viewing the columns of the obs.
Index(['assay_ontology_term_id', 'cell_type_ontology_term_id',
       'development_stage_ontology_term_id', 'disease_ontology_term_id',
       'self_reported_ethnicity_ontology_term_id', 'organism_ontology_term_id',
       'sex_ontology_term_id', 'tissue_ontology_term_id', 'is_primary_data',
       'Class', 'Subclass', 'Supertype', 'Age at death', 'donor_id',
       'suspension_type', 'Number of UMIs', 'Genes detected',
       'Fraction mitochrondrial UMIs', 'tissue_type', 'cell_type', 'assay',
       'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity',
       'development_stage', 'observation_joinid'],
      dtype='object')

### View the unique values.
Find unique values for column:  assay_ontology_term_id
['EFO:0009922']
Categories (1, object): ['EFO:0009922']
Find unique values for column:  cell_type_ontology_term_id
['CL:4023040']
Categories (1, object): ['CL:4023040']
Find unique values for column:  development_stage_ontology_term_id
['HsapDv:0000258']
Categories (1, object): ['HsapDv:0000258']
Find unique values for column:  disease_ontology_term_id
['PATO:0000461']
Categories (1, object): ['PATO:0000461']
Find unique values for column:  self_reported_ethnicity_ontology_term_id
['unknown']
Categories (1, object): ['unknown']
Find unique values for column:  organism_ontology_term_id
['NCBITaxon:9606']
Categories (1, object): ['NCBITaxon:9606']
Find unique values for column:  sex_ontology_term_id
['PATO:0000384', 'PATO:0000383']
Categories (2, object): ['PATO:0000383', 'PATO:0000384']
Find unique values for column:  tissue_ontology_term_id
['UBERON:0001950']
Categories (1, object): ['UBERON:0001950']
Find unique values for column:  is_primary_data
[ True]
Find unique values for column:  Class
['Neuronal: Glutamatergic']
Categories (1, object): ['Neuronal: Glutamatergic']
Find unique values for column:  Subclass
['L2/3 IT']
Categories (1, object): ['L2/3 IT']
Find unique values for column:  Supertype
['L2/3 IT_1', 'L2/3 IT_6', 'L2/3 IT_5', 'L2/3 IT_10', 'L2/3 IT_12', 'L2/3 IT_13', 'L2/3 IT_3', 'L2/3 IT_7', 'L2/3 IT_8']
Categories (9, object): ['L2/3 IT_1', 'L2/3 IT_6', 'L2/3 IT_7', 'L2/3 IT_5', ..., 'L2/3 IT_10',
                         'L2/3 IT_8', 'L2/3 IT_12', 'L2/3 IT_3']
Find unique values for column:  Age at death
['18 to 40 years old', '63+ years old', '41 to 62 years old']
Categories (3, object): ['18 to 40 years old', '41 to 62 years old', '63+ years old']
Find unique values for column:  donor_id
['H17.26.003', 'H18.06.363', 'H19.06.356', 'H18.03.008', 'H18.26.404', ..., 'H17.06.003', 'H19.03.302', 'H20.06.354', 'H18.03.313', 'H18.26.403']
Length: 78
Categories (78, object): ['H15.03.002', 'H15.03.003', 'H15.03.005', 'H15.03.006', ..., 'H20.03.303',
                          'H20.06.354', 'H20.26.401', 'H20.26.406']
Find unique values for column:  suspension_type
['cell']
Categories (1, object): ['cell']
Find unique values for column:  Number of UMIs
[ 8465. 16144. 29335. ... 70487. 51919. 52548.]
Find unique values for column:  Genes detected
[3638 5462 6751 ... 1874 2254 2245]
Find unique values for column:  Fraction mitochrondrial UMIs
[1.1813349e-04 6.1950195e-05 1.0227738e-04 ... 1.3483060e-03 4.9481395e-04
 1.2224260e-03]
Find unique values for column:  tissue_type
['tissue']
Categories (1, object): ['tissue']
Find unique values for column:  cell_type
['L2/3-6 intratelencephalic projecting glutamat...]
Categories (1, object): ['L2/3-6 intratelencephalic projecting glutamat...]
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
['neocortex']
Categories (1, object): ['neocortex']
Find unique values for column:  self_reported_ethnicity
['unknown']
Categories (1, object): ['unknown']
Find unique values for column:  development_stage
['adult stage']
Categories (1, object): ['adult stage']
Find unique values for column:  observation_joinid
['YkkTir=k}f' 'cH(I9*k8gl' 'EMmxgi4NVy' ... '@#zo+B2D+w' 'z6M1ylw!D*'
 'T2!~^<~#jc']

### View the adata.var:
                 highly_variable  highly_variable_rank     means  variances  \
ENSG00000227232            False                   NaN  0.039557   0.041158   
ENSG00000177757            False                   NaN  0.002662   0.002781   
ENSG00000225880            False                   NaN  0.064196   0.067636   
ENSG00000228794            False                   NaN  0.716266   1.157470   
ENSG00000230368            False                   NaN  0.003001   0.003068   

                 variances_norm  feature_is_filtered  n_cells_by_counts  \
ENSG00000227232        0.846121                False               3026   
ENSG00000177757        0.964394                False                208   
ENSG00000225880        0.823256                False               4821   
ENSG00000228794        0.875232                False              34369   
ENSG00000230368        0.941608                False                236   

                 mean_counts  log1p_mean_counts  pct_dropout_by_counts  \
ENSG00000227232     0.011550           0.011484              96.199972   
ENSG00000177757     0.000801           0.000800              99.738795   
ENSG00000225880     0.020996           0.020778              93.945825   
ENSG00000228794     0.181944           0.167160              56.839673   
ENSG00000230368     0.000961           0.000961              99.703633   

                 total_counts  log1p_total_counts feature_name  \
ENSG00000227232    919.728638            6.825165       WASH7P   
ENSG00000177757     63.763554            4.170743       FAM87B   
ENSG00000225880   1671.901733            7.422315    LINC00115   
ENSG00000228794  14488.309570            9.581166    LINC01128   
ENSG00000230368     76.551476            4.350942       FAM41C   

                feature_reference feature_biotype feature_length  \
ENSG00000227232    NCBITaxon:9606            gene           1351   
ENSG00000177757    NCBITaxon:9606            gene           1947   
ENSG00000225880    NCBITaxon:9606            gene           3312   
ENSG00000228794    NCBITaxon:9606            gene           1627   
ENSG00000230368    NCBITaxon:9606            gene            504   

                           feature_type  
ENSG00000227232  unprocessed_pseudogene  
ENSG00000177757                  lncRNA  
ENSG00000225880                  lncRNA  
ENSG00000228794                  lncRNA  
ENSG00000230368                  lncRNA  
Size of adata.var:
(18786, 17)

### Check which layer stores the raw count:
Trying to print adata.X.A[1:25, 1:25]
[[0.         0.         0.         0.         0.         0.48211867
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.48211867 0.         1.050299   0.        ]
 [0.         0.         0.5198939  0.         0.         0.
  0.29335937 0.         0.         0.5198939  0.         0.5198939
  0.         0.         0.         0.         0.         0.
  0.         0.         0.29335937 0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.4348629  0.         0.4348629
  1.1565683  0.         0.         0.         0.         0.4348629
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.3302698  0.         0.         0.         0.5781218
  0.         0.         0.         0.3302698  0.         0.
  0.5781218  0.         0.         0.         0.         0.3302698
  0.         0.         0.3302698  0.         0.3302698  0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.61470586
  0.         0.         0.         0.         0.61470586 0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.40294772
  0.40294772 0.         0.         0.         0.         0.40294772
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.68485165 0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.68485165 0.         0.         0.3999424  0.         0.68485165
  0.         0.         0.         0.         0.9063278  0.3999424 ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.20873262
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.77037555 0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.40401208
  0.40401208 0.222272   0.         0.222272   0.         0.222272
  0.         0.         0.         0.         0.69096684 0.        ]
 [0.         0.         0.         0.         0.         0.
  0.20534302 0.         0.         0.         0.         0.5210761
  0.20534302 0.         0.         0.         0.         0.20534302
  0.         0.         0.20534302 0.         0.64803725 0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.63737226 0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.         0.9085021  0.         0.         0.
  0.         0.         0.         0.         0.40114558 0.40114558
  0.         0.         0.         0.         0.         0.
  0.         0.         0.40114558 0.         0.6866609  0.        ]
 [0.         0.         0.38799572 0.         0.         0.
  0.         0.         0.         0.38799572 0.         0.
  0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.38799572 0.        ]
 [0.         0.         0.40086854 0.         0.         0.40086854
  0.         0.         0.         0.         0.         0.40086854
  0.40086854 0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.40086854 0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.4824981
  0.80657995 0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.4824981  0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.6078453  0.         0.
  0.6078453  0.         0.         0.         0.         0.34941316
  0.         0.         0.         0.         0.         0.        ]
 [0.         0.22008134 0.22008134 0.         0.         0.22008134
  0.         0.         0.         0.         0.         0.
  0.22008134 0.         0.         0.         0.         0.
  0.         0.         0.4003562  0.22008134 0.22008134 0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.36147508 0.         0.         0.         0.         0.36147508
  0.         0.         0.         0.36147508 0.62641454 0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.42581785
  0.42581785 0.         0.         0.         0.         0.
  0.         0.         0.         0.         1.2958797  0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.
  0.50706744 0.         0.         0.         0.         0.50706744
  0.         0.         0.         0.         0.9749508  0.        ]
 [0.         0.         0.5933696  0.         0.         0.
  0.         0.         0.         0.         0.34006482 0.34006482
  0.5933696  0.         0.         0.         0.         0.
  0.         0.         0.34006482 0.         0.         0.        ]
 [0.         0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.         0.7471852
  0.44181293 0.         0.         0.         0.         0.
  0.         0.         0.         0.         0.7471852  0.        ]
 [0.         0.         0.27842528 0.         0.         0.
  0.         0.         0.         0.27842528 0.27842528 0.27842528
  0.         0.         0.         0.         0.         0.27842528
  0.         0.         0.27842528 0.         0.27842528 0.        ]]
Trying to print adata.X[1:25, 1:25]
  (0, 5)	0.48211867
  (0, 20)	0.48211867
  (0, 22)	1.050299
  (1, 2)	0.5198939
  (1, 6)	0.29335937
  (1, 9)	0.5198939
  (1, 11)	0.5198939
  (1, 20)	0.29335937
  (3, 9)	0.4348629
  (3, 11)	0.4348629
  (3, 12)	1.1565683
  (3, 17)	0.4348629
  (4, 1)	0.3302698
  (4, 5)	0.5781218
  (4, 9)	0.3302698
  (4, 12)	0.5781218
  (4, 17)	0.3302698
  (4, 20)	0.3302698
  (4, 22)	0.3302698
  (5, 17)	0.61470586
  (5, 22)	0.61470586
  (6, 11)	0.40294772
  (6, 12)	0.40294772
  (6, 17)	0.40294772
  (7, 2)	0.68485165
  :	:
  (18, 12)	0.36147508
  (18, 17)	0.36147508
  (18, 21)	0.36147508
  (18, 22)	0.62641454
  (19, 11)	0.42581785
  (19, 12)	0.42581785
  (19, 22)	1.2958797
  (20, 12)	0.50706744
  (20, 17)	0.50706744
  (20, 22)	0.9749508
  (21, 2)	0.5933696
  (21, 10)	0.34006482
  (21, 11)	0.34006482
  (21, 12)	0.5933696
  (21, 20)	0.34006482
  (22, 11)	0.7471852
  (22, 12)	0.44181293
  (22, 22)	0.7471852
  (23, 2)	0.27842528
  (23, 9)	0.27842528
  (23, 10)	0.27842528
  (23, 11)	0.27842528
  (23, 17)	0.27842528
  (23, 20)	0.27842528
  (23, 22)	0.27842528
Trying to print adata.raw.X.A[1:25, 1:25]
[[0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 3. 0.]
 [0. 0. 2. 0. 0. 0. 1. 0. 0. 2. 0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 1. 4. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
 [0. 1. 0. 0. 0. 2. 0. 0. 0. 1. 0. 0. 2. 0. 0. 0. 0. 1. 0. 0. 1. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
 [0. 0. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 1. 0. 2. 0. 0. 0. 0. 3. 1.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 5. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 2. 1. 0. 1. 0. 1. 0. 0. 0. 0. 4. 0.]
 [0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 3. 1. 0. 0. 0. 0. 1. 0. 0. 1. 0. 4. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
 [0. 0. 3. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 2. 0.]
 [0. 0. 1. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 0. 1. 0. 0. 1. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 2. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 2. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
 [0. 1. 1. 0. 0. 1. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 0. 0. 0. 2. 1. 1. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0. 0. 1. 0. 0. 0. 1. 2. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 1. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 5. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0. 0. 0. 0. 2. 0. 0. 0. 0. 5. 0.]
 [0. 0. 2. 0. 0. 0. 0. 0. 0. 0. 1. 1. 2. 0. 0. 0. 0. 0. 0. 0. 1. 0. 0. 0.]
 [0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 1. 0. 0. 0. 0. 0. 0. 0. 0. 0. 2. 0.]
 [0. 0. 1. 0. 0. 0. 0. 0. 0. 1. 1. 1. 0. 0. 0. 0. 0. 1. 0. 0. 1. 0. 1. 0.]]

### View the adata.obs:
                                              assay_ontology_term_id  \
exp_component_name                                                     
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075            EFO:0009922   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075            EFO:0009922   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075            EFO:0009922   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075            EFO:0009922   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075            EFO:0009922   

                                              cell_type_ontology_term_id  \
exp_component_name                                                         
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075                 CL:4023040   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075                 CL:4023040   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075                 CL:4023040   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075                 CL:4023040   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075                 CL:4023040   

                                              development_stage_ontology_term_id  \
exp_component_name                                                                 
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075                     HsapDv:0000258   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075                     HsapDv:0000258   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075                     HsapDv:0000258   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075                     HsapDv:0000258   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075                     HsapDv:0000258   

                                              disease_ontology_term_id  \
exp_component_name                                                       
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075             PATO:0000461   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075             PATO:0000461   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075             PATO:0000461   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075             PATO:0000461   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075             PATO:0000461   

                                              self_reported_ethnicity_ontology_term_id  \
exp_component_name                                                                       
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075                                  unknown   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075                                  unknown   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075                                  unknown   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075                                  unknown   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075                                  unknown   

                                              organism_ontology_term_id  \
exp_component_name                                                        
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075            NCBITaxon:9606   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075            NCBITaxon:9606   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075            NCBITaxon:9606   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075            NCBITaxon:9606   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075            NCBITaxon:9606   

                                              sex_ontology_term_id  \
exp_component_name                                                   
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075         PATO:0000384   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075         PATO:0000384   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075         PATO:0000384   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075         PATO:0000384   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075         PATO:0000384   

                                              tissue_ontology_term_id  \
exp_component_name                                                      
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075          UBERON:0001950   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075          UBERON:0001950   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075          UBERON:0001950   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075          UBERON:0001950   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075          UBERON:0001950   

                                               is_primary_data  \
exp_component_name                                               
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075             True   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075             True   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075             True   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075             True   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075             True   

                                                                 Class  \
exp_component_name                                                       
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075  Neuronal: Glutamatergic   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075  Neuronal: Glutamatergic   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075  Neuronal: Glutamatergic   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075  Neuronal: Glutamatergic   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075  Neuronal: Glutamatergic   

                                              Subclass  Supertype  \
exp_component_name                                                  
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075  L2/3 IT  L2/3 IT_1   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075  L2/3 IT  L2/3 IT_6   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075  L2/3 IT  L2/3 IT_5   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075  L2/3 IT  L2/3 IT_1   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075  L2/3 IT  L2/3 IT_6   

                                                     Age at death    donor_id  \
exp_component_name                                                              
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075  18 to 40 years old  H17.26.003   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075  18 to 40 years old  H17.26.003   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075  18 to 40 years old  H17.26.003   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075  18 to 40 years old  H17.26.003   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075  18 to 40 years old  H17.26.003   

                                              suspension_type  Number of UMIs  \
exp_component_name                                                              
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075            cell          8465.0   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075            cell         16144.0   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075            cell         29335.0   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075            cell         11841.0   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075            cell         18359.0   

                                               Genes detected  \
exp_component_name                                              
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075            3638   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075            5462   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075            6751   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075            4345   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075            5378   

                                               Fraction mitochrondrial UMIs  \
exp_component_name                                                            
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075                      0.000118   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075                      0.000062   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075                      0.000102   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075                      0.000338   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075                      0.000109   

                                              tissue_type  \
exp_component_name                                          
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075      tissue   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075      tissue   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075      tissue   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075      tissue   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075      tissue   

                                                                                       cell_type  \
exp_component_name                                                                                 
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075  L2/3-6 intratelencephalic projecting glutamate...   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075  L2/3-6 intratelencephalic projecting glutamate...   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075  L2/3-6 intratelencephalic projecting glutamate...   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075  L2/3-6 intratelencephalic projecting glutamate...   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075  L2/3-6 intratelencephalic projecting glutamate...   

                                                   assay disease  \
exp_component_name                                                 
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075  10x 3' v3  normal   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075  10x 3' v3  normal   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075  10x 3' v3  normal   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075  10x 3' v3  normal   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075  10x 3' v3  normal   

                                                   organism   sex     tissue  \
exp_component_name                                                             
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075  Homo sapiens  male  neocortex   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075  Homo sapiens  male  neocortex   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075  Homo sapiens  male  neocortex   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075  Homo sapiens  male  neocortex   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075  Homo sapiens  male  neocortex   

                                              self_reported_ethnicity  \
exp_component_name                                                      
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075                 unknown   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075                 unknown   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075                 unknown   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075                 unknown   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075                 unknown   

                                              development_stage  \
exp_component_name                                                
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075       adult stage   
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075       adult stage   
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075       adult stage   
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075       adult stage   
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075       adult stage   

                                              observation_joinid  
exp_component_name                                                
AAACCCAAGAGAGGTA-L8TX_190416_01_C06-871132075         YkkTir=k}f  
AAACGAAGTTGACGGA-L8TX_190416_01_C06-871132075         cH(I9*k8gl  
AAACGAATCTCAAAGC-L8TX_190416_01_C06-871132075         EMmxgi4NVy  
AAACGCTGTACTCGCG-L8TX_190416_01_C06-871132075         rkaPxcDQ;Y  
AAAGAACTCTCAATCT-L8TX_190416_01_C06-871132075         QCL61AQFDF  
Size of adata.obs:
(79631, 28)

### Number of female donors: 
32

### Number of male donors: 
46