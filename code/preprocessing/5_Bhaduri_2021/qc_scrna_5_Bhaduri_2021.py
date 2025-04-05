# This script processes the scRNAseq data from prenatal neocortex from Bhaduri et al. 2021
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2024-01-09

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import os
import single_cell_helper_functions_v3
from scipy.sparse import csc_matrix
import sys
# usage: python qc_scrna_5_Bhaduri_2021.py 276_Bhaduri2021_PrenatalNeocortex_Human_2021_14wpc_PrefrontalCortex 14wpc PrefrontalCortex
def main():
# setting up #####
    id = sys.argv[1]
    dev_age = sys.argv[2]
    tissue = sys.argv[3]
    base_dir = "Preprocessing_scRNA/data" #change here if necessary
    h5ad_path = os.path.join("Preprocessing_scRNA/data/Bhaduri2021_SecondTrimester_Neocortex", tissue + ".h5ad") #change here per data
    gene_names_org = "symbol" #change here per data

    # create dev_age dict
    dev_age_dict = {"14wpc": "14th week post-fertilization human stage",
                    "16wpc": "16th week post-fertilization human stage",
                    "17wpc": "17th week post-fertilization human stage",
                    "18wpc": "18th week post-fertilization human stage",
                    "19wpc": "19th week post-fertilization human stage",
                    "20wpc": "20th week post-fertilization human stage",
                    "22wpc": "22nd week post-fertilization human stage",
                    "25wpc": "25th week post-fertilization human stage"}

    # read in adata #####
    adata = anndata.read(h5ad_path)
    metadata = adata.obs
    metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
    metadata = metadata.rename(columns = {'index': 'cell_id'})

    # since the adata is for all of the developmental age together, we will need to separate them
    print("Subsetting the adata for dev_age ", dev_age)
    log_dev_age = open(os.path.join(base_dir, id, "log.txt"), "w")
    clean_adata_dev_age_fp = os.path.join(base_dir, id, id + ".h5ad")

    # format table 2
    table2_dev_age = open(os.path.join(base_dir, id, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2_dev_age)

    # format table 3
    table3_dev_age = open(os.path.join(base_dir, id, "table3.csv"), "w")
    
    cell_id = metadata[(metadata["development_stage"] == dev_age_dict[dev_age])]["cell_id"]
    adata.obs.set_index('index', inplace=True)
    
    adata_dev_age_subset = adata[cell_id]

    print("Subsetting adata for group is finished.", file=log_dev_age)
    print("Viewing the adata observations.", file=log_dev_age)
    print(adata_dev_age_subset.obs, file=log_dev_age)

    print("Viewing the adata variables.", file=log_dev_age)
    print(adata_dev_age_subset.var, file=log_dev_age)

    print("Viewing the adata matrix - are these integer counts?",file=log_dev_age)
    print(adata_dev_age_subset.X.A[1:25,1:25],file=log_dev_age)

    adata_nrow = adata_dev_age_subset.shape[0]
    adata_ncol = adata_dev_age_subset.shape[1]
    print(f"After subsetting, adata has {adata_nrow} cells and {adata_ncol} genes.", file=log_dev_age)
    original = ["Original", str(adata_dev_age_subset.n_obs), str(adata_dev_age_subset.n_vars)]
    print(",".join(original), file=table2_dev_age)

    # mitochondrial genes
    adata_dev_age_subset.var['mt'] = adata_dev_age_subset.var["feature_name"].str.startswith("MT-")
    # count the number of genes that are mitochondrial genes
    n_mt_genes = np.count_nonzero(adata_dev_age_subset.var['mt'])
    print(f"adata has {n_mt_genes} mitochondrial genes.", file=log_dev_age)

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata_dev_age_subset, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # first filtering: accept cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    sc.pp.filter_cells(adata_dev_age_subset, min_genes=200)
    sc.pp.filter_genes(adata_dev_age_subset, min_cells=3)

    print(f"After first filtering (see documentation for definition), adata has {adata_dev_age_subset.n_obs} cells and {adata_dev_age_subset.n_vars} genes.", file=log_dev_age)
    # save to table 2
    first_filter = ["Filter#1", str(adata_dev_age_subset.n_obs), str(adata_dev_age_subset.n_vars)]
    print(",".join(first_filter), file=table2_dev_age)

    # second filtering based on mt percentage
    adata_dev_age_subset = adata_dev_age_subset[adata_dev_age_subset.obs['pct_counts_mt'] < 10, :]
    print(f"After second filtering (see documentation for definition), adata has {adata_dev_age_subset.n_obs} cells and {adata_dev_age_subset.n_vars} genes.", file=log_dev_age)
    # save to table 2
    second_filter = ["Filter#2", str(adata_dev_age_subset.n_obs), str(adata_dev_age_subset.n_vars)]
    print(",".join(second_filter), file=table2_dev_age)

    # rename cell type columns
    adata_dev_age_subset.obs.rename(columns={"cell_type": "cell_type_level_1"}, inplace = True)
    adata_dev_age_subset.obs.rename(columns={"cluster_label": "cell_type_level_2"}, inplace = True)

    # modify the adata.var so that we have the row index of adata.var is symbol and ensemble column is named ensembl
    adata_dev_age_subset.var.rename(columns={"ensg_id": "ensg_id_orig"}, inplace = True) #the row name ensg_id exists both as row name and column name so I'm renaming the column here
    adata_dev_age_subset.var.reset_index(inplace=True) #since row index of adata.var is ensemble_id, we want to convert this into a column
    adata_dev_age_subset.var = adata_dev_age_subset.var.rename(columns={'feature_name': 'symbol'}) #update the column names to have consistent naming convention
    # convert to ensembl from symbol
    gene_names_fp = 'gene_names_human.txt' #path to gene_names_human.txt
    gene_names_add = "Ensembl gene ID"
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata_dev_age_subset, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log_dev_age)
    adata_dev_age_subset.var.set_index('symbol', inplace=True) # set the index to gene symbol
    # save to table 2
    ensembl_convertable = ["Ngenes with ensembl", str(adata_dev_age_subset.n_obs), str(ngenes_after_conversion)]
    print(",".join(ensembl_convertable), file=table2_dev_age)

    # save
    adata_dev_age_subset.write_h5ad(filename=clean_adata_dev_age_fp)

    # print out the saved adata
    print(adata_dev_age_subset, file=log_dev_age)
    print(adata_dev_age_subset.var, file=log_dev_age) 
    print(adata_dev_age_subset.obs, file=log_dev_age)       

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset
    level1_cts = adata_dev_age_subset.obs["cell_type_level_1"].dropna().unique()
    for ct in level1_cts:
        ct_data = adata_dev_age_subset[adata_dev_age_subset.obs["cell_type_level_1"] == ct, :].X
        out = ["level_1", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_dev_age)
    
    level2_cts = adata_dev_age_subset.obs["cell_type_level_2"].dropna().unique()
    for ct in level2_cts:
        ct_data = adata_dev_age_subset[adata_dev_age_subset.obs["cell_type_level_2"] == ct, :].X
        out = ["level_2", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_dev_age)

    table2_dev_age.close()
    table3_dev_age.close()
    log_dev_age.close()

main()


