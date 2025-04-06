# This script processes the scRNAseq data from adult neocortex from Jorstad et al. 2023 smartseq
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2024-01-11

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import os
import single_cell_helper_functions_v3
from scipy.sparse import csc_matrix
import sys
# usage: python qc_scrna_16_Jorstad_2023_smartseq.py 326_Jorstad2023_M1_Human_2023_Smartseq M1
def main():
# setting up #####
    id = sys.argv[1]
    tissue_abbrev = sys.argv[2]
    base_dir = "Preprocessing_scRNA/data" #change here if necessary
    h5ad_path = "Jorstadetal2023_Neocortex_smartseq.h5ad" #change here per data
    gene_names_org = "symbol" #change here per data

    # create tissue dict
    tissue_dict = {"M1": "primary motor cortex", 
                   "S1": "primary somatosensory cortex",
                   "A1": "primary auditory cortex", 
                   "V1": "primary visual cortex",
                   "ACC": "anterior cingulate gyrus",
                   "MTG": "middle temporal gyrus"}

    # read in adata #####
    adata = anndata.read(h5ad_path)
    metadata = adata.obs
    metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
    metadata = metadata.rename(columns = {'index': 'cell_id'})

    # since the adata is for all of the tissue areas together, we will need to separate them
    print("Subsetting the adata for tissue", tissue_abbrev)
    log = open(os.path.join(base_dir, id, "log.txt"), "w")
    clean_adata_fp = os.path.join(base_dir, id, id + ".h5ad")

    # format table 2
    table2 = open(os.path.join(base_dir, id, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2)

    # format table 3
    table3 = open(os.path.join(base_dir, id, "table3.csv"), "w")
    
    cell_id = metadata[(metadata["tissue"] == tissue_dict[tissue_abbrev])]["cell_id"]
    adata.obs.set_index('index', inplace=True)
    
    adata_tissue_subset = adata[cell_id]

    print("Subsetting adata for group is finished.", file=log)
    print("Viewing the adata observations.", file=log)
    print(adata_tissue_subset.obs, file=log)

    print("Viewing the adata variables.", file=log)
    print(adata_tissue_subset.var, file=log)

    print("Viewing the adata matrix - are these integer counts?",file=log)
    print(adata_tissue_subset.X.A[1:25,1:25],file=log)

    adata_nrow = adata_tissue_subset.shape[0]
    adata_ncol = adata_tissue_subset.shape[1]
    print(f"After subsetting, adata has {adata_nrow} cells and {adata_ncol} genes.", file=log)
    original = ["Original", str(adata_tissue_subset.n_obs), str(adata_tissue_subset.n_vars)]
    print(",".join(original), file=table2)

    # mitochondrial genes
    adata_tissue_subset.var['mt'] = adata_tissue_subset.var["feature_name"].str.startswith("MT-")
    # count the number of genes that are mitochondrial genes
    n_mt_genes = np.count_nonzero(adata_tissue_subset.var['mt'])
    print(f"adata has {n_mt_genes} mitochondrial genes.", file=log)

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata_tissue_subset, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # first filtering: accept cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    sc.pp.filter_cells(adata_tissue_subset, min_genes=200)
    sc.pp.filter_genes(adata_tissue_subset, min_cells=3)

    print(f"After first filtering (see documentation for definition), adata has {adata_tissue_subset.n_obs} cells and {adata_tissue_subset.n_vars} genes.", file=log)
    # save to table 2
    first_filter = ["Filter#1", str(adata_tissue_subset.n_obs), str(adata_tissue_subset.n_vars)]
    print(",".join(first_filter), file=table2)

    # second filtering based on mt percentage
    adata_tissue_subset = adata_tissue_subset[adata_tissue_subset.obs['pct_counts_mt'] < 10, :]
    print(f"After second filtering (see documentation for definition), adata has {adata_tissue_subset.n_obs} cells and {adata_tissue_subset.n_vars} genes.", file=log)
    # save to table 2
    second_filter = ["Filter#2", str(adata_tissue_subset.n_obs), str(adata_tissue_subset.n_vars)]
    print(",".join(second_filter), file=table2)

    # rename cell type columns
    adata_tissue_subset.obs.rename(columns={"class": "cell_type_level_1"}, inplace = True)
    adata_tissue_subset.obs.rename(columns={"subclass": "cell_type_level_2"}, inplace = True)
    adata_tissue_subset.obs.rename(columns={"cluster": "cell_type_level_3"}, inplace = True)

    # modify the adata.var so that we have the row index of adata.var is symbol and ensemble column is named ensembl
    adata_tissue_subset.var.rename(columns={"ensg_id": "ensg_id_orig"}, inplace = True) #the row name ensg_id exists both as row name and column name so I'm renaming the column here
    adata_tissue_subset.var.reset_index(inplace=True) #since row index of adata.var is ensemble_id, we want to convert this into a column
    adata_tissue_subset.var = adata_tissue_subset.var.rename(columns={'feature_name': 'symbol'}) #update the column names to have consistent naming convention
    # convert to ensembl from symbol
    gene_names_fp = 'gene_names_human.txt' #path to gene_names_human.txt
    gene_names_add = "Ensembl gene ID"
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata_tissue_subset, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log)
    adata_tissue_subset.var.set_index('symbol', inplace=True) # set the index to gene symbol
    # save to table 2
    ensembl_convertable = ["Ngenes with ensembl", str(adata_tissue_subset.n_obs), str(ngenes_after_conversion)]
    print(",".join(ensembl_convertable), file=table2)

    # save
    adata_tissue_subset.write_h5ad(filename=clean_adata_fp)

    # print out the saved adata
    print(adata_tissue_subset, file=log)
    print(adata_tissue_subset.var, file=log) 
    print(adata_tissue_subset.obs, file=log)       

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset
    level1_cts = adata_tissue_subset.obs["cell_type_level_1"].dropna().unique()
    for ct in level1_cts:
        ct_data = adata_tissue_subset[adata_tissue_subset.obs["cell_type_level_1"] == ct, :].X
        out = ["level_1", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)
    
    level2_cts = adata_tissue_subset.obs["cell_type_level_2"].dropna().unique()
    for ct in level2_cts:
        ct_data = adata_tissue_subset[adata_tissue_subset.obs["cell_type_level_2"] == ct, :].X
        out = ["level_2", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)
    
    level3_cts = adata_tissue_subset.obs["cell_type_level_3"].dropna().unique()
    for ct in level3_cts:
        ct_data = adata_tissue_subset[adata_tissue_subset.obs["cell_type_level_3"] == ct, :].X
        out = ["level_3", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)

    table2.close()
    table3.close()
    log.close()

main()