# This script processes the scRNAseq data from Braun et al. 2023
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2024-02-11

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import os
import single_cell_helper_functions_v3
from scipy.sparse import csc_matrix
import sys

# usage: python qc_scrna_23_Braun_2023.py 443_Braun2023_Human_2023_FirstTrimester_Brain_CarnegieStage18 CarnegieStage18 Brain

def main():
# setting up #####
    id = sys.argv[1]
    tissue = sys.argv[2]
    dev_age = sys.argv[3]
    base_dir = "Preprocessing_scRNA/data" #change here if necessary
    h5ad_path = os.path.join("Preprocessing_scRNA/data/Braun2023_FirstTrimesterBrain", tissue + "_" + dev_age + ".h5ad") #change here per data
    gene_names_org = "symbol" #change here per data

    # read in adata #####
    adata = anndata.read(h5ad_path)
    adata.obs_names_make_unique
    # make the X layer to have the raw count layer

    log = open(os.path.join(base_dir, id, "log.txt"), "w")
    clean_adata_fp = os.path.join(base_dir, id, id + ".h5ad")

    # format table 2
    table2 = open(os.path.join(base_dir, id, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2)

    # format table 3
    table3 = open(os.path.join(base_dir, id, "table3.csv"), "w")
    
    print("Viewing the adata observations.", file=log)
    print(adata.obs, file=log)

    print("Viewing the adata variables.", file=log)
    print(adata.var, file=log)

    print("Viewing the adata matrix - are these integer counts?",file=log)
    print(adata.X.A[1:25,1:25],file=log)

    adata_nrow = adata.shape[0]
    adata_ncol = adata.shape[1]
    print(f"adata has {adata_nrow} cells and {adata_ncol} genes.", file=log)
    original = ["Original", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(original), file=table2)

    # load in the var file
    var_file = pd.read_csv("Preprocessing_scRNA/data/Braun2023_FirstTrimesterBrain/adata_var.txt", sep="\t")
    adata.var = var_file
    # mitochondrial genes
    adata.var['mt'] = adata.var["feature_name"].str.startswith("MT-")
    # count the number of genes that are mitochondrial genes
    n_mt_genes = np.count_nonzero(adata.var['mt'])
    print(f"adata has {n_mt_genes} mitochondrial genes.", file=log)

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # first filtering: accept cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    print(f"After first filtering (see documentation for definition), adata has {adata.n_obs} cells and {adata.n_vars} genes.", file=log)
    # save to table 2
    first_filter = ["Filter#1", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(first_filter), file=table2)

    # second filtering based on mt percentage
    adata = adata[adata.obs['pct_counts_mt'] < 10, :]
    print(f"After second filtering (see documentation for definition), adata has {adata.n_obs} cells and {adata.n_vars} genes.", file=log)
    # save to table 2
    second_filter = ["Filter#2", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(second_filter), file=table2)

    # rename cell type columns
    adata.obs.rename(columns={"cell_type": "cell_type_level_1"}, inplace = True)


    adata.var = adata.var.rename(columns={'feature_name': 'symbol'}) #update the column names to have consistent naming convention
    # convert to ensembl from symbol
    gene_names_fp = 'gene_names_human.txt' #path to gene_names_human.txt
    gene_names_add = "Ensembl gene ID"
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log)
    adata.var.set_index('symbol', inplace=True) # set the index to gene symbol
    # save to table 2
    ensembl_convertable = ["Ngenes with ensembl", str(adata.n_obs), str(ngenes_after_conversion)]
    print(",".join(ensembl_convertable), file=table2)

    # save
    adata.write_h5ad(filename=clean_adata_fp)

    # print out the saved adata
    print(adata, file=log)
    print(adata.var, file=log) 
    print(adata.obs, file=log)       

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset
    level1_cts = adata.obs["cell_type_level_1"].dropna().unique()
    for ct in level1_cts:
        ct_data = adata[adata.obs["cell_type_level_1"] == ct, :].X
        out = ["level_1", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)

    table2.close()
    table3.close()
    log.close()

main()


