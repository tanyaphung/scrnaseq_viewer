# This script processes the scRNAseq data from neocortex from Zhu et al. 2023 group_1
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2023-12-24

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import os
import single_cell_helper_functions_v3
from scipy.sparse import csc_matrix

def main():
    # setting up #####
    base_dir = "Preprocessing_scRNA/data" #change here if necessary
    h5ad_path = "2f17c183-388a-4c08-9adb-a146833e57ab.h5ad" #change here per data
    gene_names_org = "symbol" #change here per data

    # read in adata #####
    adata = anndata.read(h5ad_path)
    # make the X layer to have the raw count layer
    print("Make the X layer to have the raw count layer")
    adata.X = adata.raw.X.copy()
    print("Finished with making the X layer to have the raw count layer")
    metadata = adata.obs
    metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
    metadata = metadata.rename(columns = {'index': 'cell_id'})

    # since the adata is for all of the ages together, we will need to separate them

    # --------------------------------------
    # processing for group_1
    # --------------------------------------
    print("Subsetting the adata for sequencing group_1")
    id_group1 = "247_Zhu2023_Neocortex_Human_2023_group1"
    log_group1 = open(os.path.join(base_dir, id_group1, "log.txt"), "w")
    clean_adata_group1_fp = os.path.join(base_dir, id_group1, id_group1 + ".h5ad")

    # format table 2
    table2_group1 = open(os.path.join(base_dir, id_group1, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2_group1)

    # format table 3
    table3_group1 = open(os.path.join(base_dir, id_group1, "table3.csv"), "w")
    
    cell_id_group1 = metadata[(metadata["age_group"] == "early fetal")]["cell_id"]
    adata.obs.set_index('index', inplace=True)
    
    adata_group1_subset = adata[cell_id_group1]

    print("Subsetting adata for group is finished.", file=log_group1)
    print("Viewing the adata observations.", file=log_group1)
    print(adata_group1_subset.obs, file=log_group1)

    print("Viewing the adata variables.", file=log_group1)
    print(adata_group1_subset.var, file=log_group1)

    print("Viewing the adata matrix - are these integer counts?",file=log_group1)
    print(adata_group1_subset.X.A[1:25,1:25],file=log_group1)

    adata_nrow = adata_group1_subset.shape[0]
    adata_ncol = adata_group1_subset.shape[1]
    print(f"After subsetting, adata has {adata_nrow} cells and {adata_ncol} genes.", file=log_group1)
    original = ["Original", str(adata_group1_subset.n_obs), str(adata_group1_subset.n_vars)]
    print(",".join(original), file=table2_group1)

    # mitochondrial genes
    adata_group1_subset.var['mt'] = adata_group1_subset.var["feature_name"].str.startswith("MT-")
    # count the number of genes that are mitochondrial genes
    n_mt_genes = np.count_nonzero(adata_group1_subset.var['mt'])
    print(f"adata has {n_mt_genes} mitochondrial genes.", file=log_group1)

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata_group1_subset, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # first filtering: accept cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    sc.pp.filter_cells(adata_group1_subset, min_genes=200)
    sc.pp.filter_genes(adata_group1_subset, min_cells=3)

    print(f"After first filtering (see documentation for definition), adata has {adata_group1_subset.n_obs} cells and {adata_group1_subset.n_vars} genes.", file=log_group1)
    # save to table 2
    first_filter = ["Filter#1", str(adata_group1_subset.n_obs), str(adata_group1_subset.n_vars)]
    print(",".join(first_filter), file=table2_group1)

    # second filtering based on mt percentage
    adata_group1_subset = adata_group1_subset[adata_group1_subset.obs['pct_counts_mt'] < 10, :]
    print(f"After second filtering (see documentation for definition), adata has {adata_group1_subset.n_obs} cells and {adata_group1_subset.n_vars} genes.", file=log_group1)
    # save to table 2
    second_filter = ["Filter#2", str(adata_group1_subset.n_obs), str(adata_group1_subset.n_vars)]
    print(",".join(second_filter), file=table2_group1)

    # rename cell type columns
    adata_group1_subset.obs.rename(columns={"cell_type": "cell_type_level_1"}, inplace = True)
    adata_group1_subset.obs.rename(columns={"author_cell_type": "cell_type_level_2"}, inplace = True)

    # modify the adata.var so that we have the row index of adata.var is symbol and ensemble column is named ensembl
    adata_group1_subset.var.reset_index(inplace=True) #since row index of adata.var is ensemble_id_group1s, we want to convert this into a column
    adata_group1_subset.var = adata_group1_subset.var.rename(columns={'feature_name': 'symbol'}) #update the column names to have consistent naming convention
    # convert to ensembl from symbol
    gene_names_fp = 'gene_names_human.txt' #path to gene_names_human.txt
    gene_names_add = "Ensembl gene ID"
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata_group1_subset, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log_group1)
    adata_group1_subset.var.set_index('symbol', inplace=True) # set the index to gene symbol
    # save to table 2
    ensembl_convertable = ["Ngenes with ensembl", str(adata_group1_subset.n_obs), str(ngenes_after_conversion)]
    print(",".join(ensembl_convertable), file=table2_group1)

    # save
    adata_group1_subset.write_h5ad(filename=clean_adata_group1_fp)

    # print out the saved adata
    print(adata_group1_subset, file=log_group1)
    print(adata_group1_subset.var, file=log_group1) 
    print(adata_group1_subset.obs, file=log_group1)       

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset
    level1_cts = adata_group1_subset.obs["cell_type_level_1"].dropna().unique()
    for ct in level1_cts:
        ct_data = adata_group1_subset[adata_group1_subset.obs["cell_type_level_1"] == ct, :].X
        out = ["level_1", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_group1)
    level2_cts = adata_group1_subset.obs["cell_type_level_2"].dropna().unique()
    for ct in level2_cts:
        ct_data = adata_group1_subset[adata_group1_subset.obs["cell_type_level_2"] == ct, :].X
        out = ["level_2", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_group1)

    table2_group1.close()
    table3_group1.close()
    log_group1.close()

main()


