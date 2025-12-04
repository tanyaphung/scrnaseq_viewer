# Documentation
## General information
- Link: https://cellxgene.cziscience.com/collections/7d66d871-091f-4602-9f42-85f86d2853e0
- https://github.com/linnarsson-lab/human-meninges-development

- Raw counts downloaded on 2025-09-03

```
wget https://datasets.cellxgene.cziscience.com/f82c6e2b-8056-45a3-967a-fbadf057566f.h5ad
```

- Subset by tissues
'meninx', 'brain meninx', 'forebrain meninges'

- Note that there are 9 developmental stages ('10th week post-fertilization stage', '9th week post-fertilization stage', '13th week post-fertilization stage', '12th week post-fertilization stage', 'Carnegie stage 20', 'Carnegie stage 18', 'Carnegie stage 21', 'Carnegie stage 16', 'Carnegie stage 23'). I'm not separating this and just keep as prenatal. 

```
python qc_scrna_Linnarsson_meninx.py 554_Linnarsson_Prenatal_Meninx meninx
```

```
python qc_scrna_Linnarsson_meninx.py 555_Linnarsson_Prenatal_BrainMeninx brainmeninx
```

```
python qc_scrna_Linnarsson_meninx.py 556_Linnarsson_Prenatal_ForebrainMeninges forebrainmeninges
```
