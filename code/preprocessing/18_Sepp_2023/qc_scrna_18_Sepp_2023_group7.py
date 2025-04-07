# This script processes the scRNAseq data from cerebellum from Sepp et al. 2023 Nature group_7
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2023-12-22

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
    h5ad_path = "6b141abe-83b3-439d-a66d-a3b5495c52c7.h5ad" #change here per data
    gene_names_org = "symbol" #change here per data

    # read in adata #####
    adata = anndata.read(h5ad_path)

    metadata = adata.obs
    metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
    metadata = metadata.rename(columns = {'index': 'cell_id'})

    # since the adata is for all of the ages together, we will need to separate them

    # --------------------------------------
    # processing for group_7
    # --------------------------------------
    print("Subsetting the adata for sequencing group_7")
    id_group7 = "242_Sepp2023_Cerebellum_Human_2023_group7"
    log_group7 = open(os.path.join(base_dir, id_group7, "log.txt"), "w")
    clean_adata_group7_fp = os.path.join(base_dir, id_group7, id_group7 + ".h5ad")

    # format table 2
    table2_group7 = open(os.path.join(base_dir, id_group7, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2_group7)

    # format table 3
    table3_group7 = open(os.path.join(base_dir, id_group7, "table3.csv"), "w")

    normal_id = metadata[(metadata["disease"] == "normal")]["cell_id"]
    adata.obs.set_index('index', inplace=True)

    # subset for normal
    adata_normal = adata[normal_id]
    metadata_normal = adata_normal.obs
    metadata_normal.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
    metadata_normal = metadata_normal.rename(columns = {'index': 'cell_id'})
    cell_id_group7 = metadata_normal[(metadata_normal["author_stage"] == "20 wpc")]["cell_id"]
    adata_normal.obs.set_index('index', inplace=True)
    adata_group7_subset = adata_normal[cell_id_group7]

    print("Subsetting adata for group is finished.", file=log_group7)
    print("Viewing the adata observations.", file=log_group7)
    print(adata_group7_subset.obs, file=log_group7)

    print("Viewing the adata variables.", file=log_group7)
    print(adata_group7_subset.var, file=log_group7)

    print("Viewing the adata matrix - are these integer counts?",file=log_group7)
    print(adata_group7_subset.X.A[1:10,1:10],file=log_group7)

    adata_nrow = adata_group7_subset.shape[0]
    adata_ncol = adata_group7_subset.shape[1]
    print(f"After subsetting, adata has {adata_nrow} cells and {adata_ncol} genes.", file=log_group7)
    original = ["Original", str(adata_group7_subset.n_obs), str(adata_group7_subset.n_vars)]
    print(",".join(original), file=table2_group7)

    # mitochondrial genes
    adata_group7_subset.var['mt'] = adata_group7_subset.var["feature_name"].str.startswith("MT-")
    # count the number of genes that are mitochondrial genes
    n_mt_genes = np.count_nonzero(adata_group7_subset.var['mt'])
    print(f"adata has {n_mt_genes} mitochondrial genes.", file=log_group7)

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata_group7_subset, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # first filtering: accept cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    sc.pp.filter_cells(adata_group7_subset, min_genes=200)
    sc.pp.filter_genes(adata_group7_subset, min_cells=3)

    print(f"After first filtering (see documentation for definition), adata has {adata_group7_subset.n_obs} cells and {adata_group7_subset.n_vars} genes.", file=log_group7)
    # save to table 2
    first_filter = ["Filter#1", str(adata_group7_subset.n_obs), str(adata_group7_subset.n_vars)]
    print(",".join(first_filter), file=table2_group7)

    # second filtering based on mt percentage
    adata_group7_subset = adata_group7_subset[adata_group7_subset.obs['pct_counts_mt'] < 10, :]
    print(f"After second filtering (see documentation for definition), adata has {adata_group7_subset.n_obs} cells and {adata_group7_subset.n_vars} genes.", file=log_group7)
    # save to table 2
    second_filter = ["Filter#2", str(adata_group7_subset.n_obs), str(adata_group7_subset.n_vars)]
    print(",".join(second_filter), file=table2_group7)

    # rename cell type columns
    adata_group7_subset.obs.rename(columns={"cell_type": "cell_type_level_1"}, inplace = True)

    # modify the adata.var so that we have the row index of adata.var is symbol and ensemble column is named ensembl
    adata_group7_subset.var.reset_index(inplace=True) #since row index of adata.var is ensemble_id_group7s, we want to convert this into a column
    adata_group7_subset.var = adata_group7_subset.var.rename(columns={'feature_name': 'symbol'}) #update the column names to have consistent naming convention
    # convert to ensembl from symbol
    gene_names_fp = 'gene_names_human.txt' #path to gene_names_human.txt
    gene_names_add = "Ensembl gene ID"
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata_group7_subset, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log_group7)
    adata_group7_subset.var.set_index('symbol', inplace=True) # set the index to gene symbol
    # save to table 2
    ensembl_convertable = ["Ngenes with ensembl", str(adata_group7_subset.n_obs), str(ngenes_after_conversion)]
    print(",".join(ensembl_convertable), file=table2_group7)

    # save
    adata_group7_subset.write_h5ad(filename=clean_adata_group7_fp)

    # print out the saved adata
    print(adata_group7_subset, file=log_group7)
    print(adata_group7_subset.var, file=log_group7) 
    print(adata_group7_subset.obs, file=log_group7)       

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset
    level1_cts = adata_group7_subset.obs["cell_type_level_1"].dropna().unique()
    for ct in level1_cts:
        ct_data = adata_group7_subset[adata_group7_subset.obs["cell_type_level_1"] == ct, :].X
        out = ["level_1", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_group7)

    table2_group7.close()
    table3_group7.close()
    log_group7.close()

main()


