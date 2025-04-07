# Rationale: 
- In Xu et al., they developed a methods called CellHint to harmonized scRNAseq data. Therefore, data from different studies are harmonized and cell type annotation was performed on the harmonized data. 
- In the column `Dataset` from `adata.obs`, it lists all the studies that were used for scRNAseq harmonization. 
- Typically, I would separate out the different studies. However, because in this study scRNAseq data from different studies were harmonized. Therefore, for the datasets from Xu et al., I am not separating out by datasets. 
  - If the tissue has multiple regions, I separate out by regions. For example, there are 6 regions for Heart. 

# Documentation for Blood
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/369e3ca7-5e0f-417a-ac06-641e9f274b10.h5ad
mv 369e3ca7-5e0f-417a-ac06-641e9f274b10.h5ad Blood.h5ad
```

# Documentation for Bone_marrow
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/b2eca8f3-b461-45fd-8639-890bbbf050aa.h5ad
mv b2eca8f3-b461-45fd-8639-890bbbf050aa.h5ad Bone_marrow.h5ad
```

# Documentation for Heart
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/7f37f4c6-6c72-4896-9f79-edd901ffee90.h5ad
mv 7f37f4c6-6c72-4896-9f79-edd901ffee90.h5ad Heart.h5ad
```

# Documentation for Hippocampus
- Link: https://cellxgene.cziscience.com/collections/854c0855-23ad-4362-8b77-6b1639e7a9fc
- Raw counts downloaded on 2025-03-12

```
wget https://datasets.cellxgene.cziscience.com/3be1614f-e391-47ea-9bc1-0c65eba0c540.h5ad
mv 3be1614f-e391-47ea-9bc1-0c65eba0c540.h5ad Hippocampus.h5ad
```

# Documentation for Intestine
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/4906f3d2-16f5-4498-86aa-78971994eba1.h5ad
mv 4906f3d2-16f5-4498-86aa-78971994eba1.h5ad Intestine.h5ad
```

# Documentation for Kidney
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/9f570ea4-9de1-4de4-9bfc-e6a80ad8184c.h5ad
mv 9f570ea4-9de1-4de4-9bfc-e6a80ad8184c.h5ad Kidney.h5ad
```

# Documentation for Liver
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/03107eab-699e-4c39-977c-25d685345a33.h5ad
mv 03107eab-699e-4c39-977c-25d685345a33.h5ad Liver.h5ad
```

# Documentation for Lung
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/5311ca08-a915-4bea-a83c-5f2231ba18ef.h5ad
mv 5311ca08-a915-4bea-a83c-5f2231ba18ef.h5ad Lung.h5ad
```

# Documentation for Lymph_node
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/605ef1ae-3845-44cc-b056-266eebd96232.h5ad
mv 605ef1ae-3845-44cc-b056-266eebd96232.h5ad Lymph_node.h5ad
```

# Documentation for Pancreas
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/c2fe4e46-ca49-44cc-9237-eb73b7a0b36a.h5ad
mv c2fe4e46-ca49-44cc-9237-eb73b7a0b36a.h5ad Pancreas.h5ad
```

# Documentation for Skeletal_muscle
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/9c127049-9eb8-4d2b-bb2d-ef8714c583dd.h5ad
mv 9c127049-9eb8-4d2b-bb2d-ef8714c583dd.h5ad Skeletal_muscle.h5ad
```

# Documentation for Spleen
- Raw counts downloaded on 2025-04-01

```
wget https://datasets.cellxgene.cziscience.com/80012b95-462d-4e27-abdc-12a64b1081f3.h5ad
mv 80012b95-462d-4e27-abdc-12a64b1081f3.h5ad Spleen.h5ad
```







