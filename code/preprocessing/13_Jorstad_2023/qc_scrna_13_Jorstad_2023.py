# This script processes the scRNAseq data for MTG from Jorstad et al. 2023 Science
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2023-12-18

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
# from matplotlib import pyplot as plt
import os
import single_cell_helper_functions_v3
from scipy.sparse import csc_matrix

def main():
    # setting up #####
    # ct_colname = "Subclass" #change here per data
    base_dir = "Preprocessing_scRNA/data" #change here if necessary
    h5ad_path = "fafef79b-a863-4e63-a162-dd933d552ace.h5ad" #change here per data
    gene_names_org = "symbol" #change here per data

    # read in adata #####
    adata = anndata.read(h5ad_path)

    # make the X layer to have the raw count layer
    print("Make the X layer to have the raw count layer")
    adata.X = adata.raw.X.copy()
    print("Finished with making the X layer to have the raw count layer")

    # adata.obs.set_index('index', inplace=True)
    metadata = adata.obs
    metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
    metadata = metadata.rename(columns = {'index': 'cell_id'})

    # since the adata is for all of the sequencing methods together, we will need to separate them #####
    # processing for smartseq
    print("Subsetting the adata for sequencing method #1 Smart-seq v4")
    id_smartseq = "226_Jorstadetal2023_MTG_Human_2023_smartseq"
    log_smartseq = open(os.path.join(base_dir, id_smartseq, "log.txt"), "w")
    clean_adata_smartseq_fp = os.path.join(base_dir, id_smartseq, id_smartseq + ".h5ad")

    # format table 2
    table2_smartseq = open(os.path.join(base_dir, id_smartseq, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2_smartseq)

    # format table 3
    table3_smartseq = open(os.path.join(base_dir, id_smartseq, "table3.csv"), "w")
    
    cell_id_smartseqs = metadata[metadata["assay"] == "Smart-seq v4"]["cell_id"]
    adata.obs.set_index('index', inplace=True)
    
    adata_smartseq_subset = adata[cell_id_smartseqs]
    # adata_smartseq_subset.X = adata_smartseq_subset.raw.X.copy() #copy the raw layer to X layer

    print("Subsetting adata for smart seq is finished.", file=log_smartseq)
    print("Viewing the adata observations.", file=log_smartseq)
    print(adata_smartseq_subset.obs, file=log_smartseq)

    print("Viewing the adata variables.", file=log_smartseq)
    print(adata_smartseq_subset.var, file=log_smartseq)

    print("Viewing the adata matrix - are these integer counts?",file=log_smartseq)
    print(adata_smartseq_subset.X[1:10,1:10],file=log_smartseq)
    print(adata_smartseq_subset.raw.X.A[1:10,1:10],file=log_smartseq)

    adata_nrow = adata_smartseq_subset.shape[0]
    adata_ncol = adata_smartseq_subset.shape[1]
    print(f"After subsetting, adata has {adata_nrow} cells and {adata_ncol} genes.", file=log_smartseq)
    original = ["Original", str(adata_smartseq_subset.n_obs), str(adata_smartseq_subset.n_vars)]
    print(",".join(original), file=table2_smartseq)

    # mitochondrial genes
    adata_smartseq_subset.var['mt'] = adata_smartseq_subset.var["feature_name"].str.startswith("MT-")
    # count the number of genes that are mitochondrial genes
    n_mt_genes = np.count_nonzero(adata_smartseq_subset.var['mt'])
    print(f"adata has {n_mt_genes} mitochondrial genes.", file=log_smartseq)

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata_smartseq_subset, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # first filtering: accept cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    sc.pp.filter_cells(adata_smartseq_subset, min_genes=200)
    sc.pp.filter_genes(adata_smartseq_subset, min_cells=3)

    print(f"After first filtering (see documentation for definition), adata has {adata_smartseq_subset.n_obs} cells and {adata_smartseq_subset.n_vars} genes.", file=log_smartseq)
    # save to table 2
    first_filter = ["Filter#1", str(adata_smartseq_subset.n_obs), str(adata_smartseq_subset.n_vars)]
    print(",".join(first_filter), file=table2_smartseq)

    # second filtering based on mt percentage
    adata_smartseq_subset = adata_smartseq_subset[adata_smartseq_subset.obs['pct_counts_mt'] < 10, :]
    print(f"After second filtering (see documentation for definition), adata has {adata_smartseq_subset.n_obs} cells and {adata_smartseq_subset.n_vars} genes.", file=log_smartseq)
    # save to table 2
    second_filter = ["Filter#2", str(adata_smartseq_subset.n_obs), str(adata_smartseq_subset.n_vars)]
    print(",".join(second_filter), file=table2_smartseq)

    # rename cell type columns
    adata_smartseq_subset.obs.rename(columns={"cell_type": "cell_type_level_1", "Subclass": "cell_type_level_2", "Cluster": "cell_type_level_3"}, inplace = True)

    # modify the adata.var so that we have the row index of adata.var is symbol and ensemble column is named ensembl
    adata_smartseq_subset.var.reset_index(inplace=True) #since row index of adata.var is ensemble_id_smartseqs, we want to convert this into a column
    adata_smartseq_subset.var = adata_smartseq_subset.var.rename(columns={'feature_name': 'symbol'}) #update the column names to have consistent naming convention
    # convert to ensembl from symbol
    gene_names_fp = 'gene_names_human.txt' #path to gene_names_human.txt
    gene_names_add = "Ensembl gene ID"
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata_smartseq_subset, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log_smartseq)
    adata_smartseq_subset.var.set_index('symbol', inplace=True) # set the index to gene symbol
    # save to table 2
    ensembl_convertable = ["Ngenes with ensembl", str(adata_smartseq_subset.n_obs), str(ngenes_after_conversion)]
    print(",".join(ensembl_convertable), file=table2_smartseq)

    # save
    adata_smartseq_subset.write_h5ad(filename=clean_adata_smartseq_fp)

    # print out the saved adata
    print(adata_smartseq_subset, file=log_smartseq)
    print(adata_smartseq_subset.var, file=log_smartseq) 
    print(adata_smartseq_subset.obs, file=log_smartseq)       

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset
    level1_cts = adata_smartseq_subset.obs["cell_type_level_1"].dropna().unique()
    for ct in level1_cts:
        ct_data = adata_smartseq_subset[adata_smartseq_subset.obs["cell_type_level_1"] == ct, :].X
        out = ["level_1", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_smartseq)

    level2_cts = adata_smartseq_subset.obs["cell_type_level_2"].dropna().unique()
    for ct in level2_cts:
        ct_data = adata_smartseq_subset[adata_smartseq_subset.obs["cell_type_level_2"] == ct, :].X
        out = ["level_2", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_smartseq)

    level3_cts = adata_smartseq_subset.obs["cell_type_level_3"].dropna().unique()
    for ct in level3_cts:
        ct_data = adata_smartseq_subset[adata_smartseq_subset.obs["cell_type_level_3"] == ct, :].X
        out = ["level_3", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_smartseq)

    table2_smartseq.close()
    table3_smartseq.close()
    log_smartseq.close()

    # processing for 10x
    # read in adata #####
    adata = anndata.read(h5ad_path)
        # make the X layer to have the raw count layer
    adata.X = adata.raw.X.copy()

    # adata.obs.set_index('index', inplace=True)
    metadata = adata.obs
    metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
    metadata = metadata.rename(columns = {'index': 'cell_id'})
    id_10x = "227_Jorstadetal2023_MTG_Human_2023_10x"
    log_10x = open(os.path.join(base_dir, id_10x, "log.txt"), "w")
    print("Subsetting the adata for sequencing method #2 10x", file=log_10x)
    clean_adata_10x_fp = os.path.join(base_dir, id_10x, id_10x + ".h5ad")

    # format table 2
    table2_10x = open(os.path.join(base_dir, id_10x, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2_10x)

    # format table 3
    table3_10x = open(os.path.join(base_dir, id_10x, "table3.csv"), "w")
    
    cell_id_10xs = metadata[metadata["assay"] == "10x 3' v3"]["cell_id"]
    adata.obs.set_index('index', inplace=True)
    
    adata_10x_subset = adata[cell_id_10xs]
    # adata_10x_subset.X.A = adata_10x_subset.raw.X.A.copy() #copy the raw layer to X layer

    print("Subsetting adata for smart seq is finished.", file=log_10x)
    print("Viewing the adata observations.", file=log_10x)
    print(adata_10x_subset.obs, file=log_10x)

    print("Viewing the adata variables.", file=log_10x)
    print(adata_10x_subset.var, file=log_10x)

    print("Viewing the adata matrix - are these integer counts?",file=log_10x)
    print(adata_10x_subset.raw.X.A[1:10,1:10],file=log_10x)

    adata_nrow = adata_10x_subset.shape[0]
    adata_ncol = adata_10x_subset.shape[1]
    print(f"After subsetting, adata has {adata_nrow} cells and {adata_ncol} genes.", file=log_10x)
    original = ["Original", str(adata_10x_subset.n_obs), str(adata_10x_subset.n_vars)]
    print(",".join(original), file=table2_10x)

    # mitochondrial genes
    adata_10x_subset.var['mt'] = adata_10x_subset.var["feature_name"].str.startswith("MT-")
    # count the number of genes that are mitochondrial genes
    n_mt_genes = np.count_nonzero(adata_10x_subset.var['mt'])
    print(f"adata has {n_mt_genes} mitochondrial genes.", file=log_10x)

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata_10x_subset, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # first filtering: accept cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    sc.pp.filter_cells(adata_10x_subset, min_genes=200)
    sc.pp.filter_genes(adata_10x_subset, min_cells=3)

    print(f"After first filtering (see documentation for definition), adata has {adata_10x_subset.n_obs} cells and {adata_10x_subset.n_vars} genes.", file=log_10x)
    # save to table 2
    first_filter = ["Filter#1", str(adata_10x_subset.n_obs), str(adata_10x_subset.n_vars)]
    print(",".join(first_filter), file=table2_10x)

    # second filtering based on mt percentage
    adata_10x_subset = adata_10x_subset[adata_10x_subset.obs['pct_counts_mt'] < 10, :]
    print(f"After second filtering (see documentation for definition), adata has {adata_10x_subset.n_obs} cells and {adata_10x_subset.n_vars} genes.", file=log_10x)
    # save to table 2
    second_filter = ["Filter#2", str(adata_10x_subset.n_obs), str(adata_10x_subset.n_vars)]
    print(",".join(second_filter), file=table2_10x)

    # rename cell type columns
    adata_10x_subset.obs.rename(columns={"cell_type": "cell_type_level_1", "Subclass": "cell_type_level_2", "Cluster": "cell_type_level_3"}, inplace = True)

    # modify the adata.var so that we have the row index of adata.var is symbol and ensemble column is named ensembl
    adata_10x_subset.var.reset_index(inplace=True) #since row index of adata.var is ensemble_id_10xs, we want to convert this into a column
    adata_10x_subset.var = adata_10x_subset.var.rename(columns={'feature_name': 'symbol'}) #update the column names to have consistent naming convention
    # convert to ensembl from symbol
    gene_names_fp = 'gene_names_human.txt' #path to gene_names_human.txt
    gene_names_add = "Ensembl gene ID"
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata_10x_subset, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log_10x)
    adata_10x_subset.var.set_index('symbol', inplace=True) # set the index to gene symbol
    # save to table 2
    ensembl_convertable = ["Ngenes with ensembl", str(adata_10x_subset.n_obs), str(ngenes_after_conversion)]
    print(",".join(ensembl_convertable), file=table2_10x)

    # save
    adata_10x_subset.write_h5ad(filename=clean_adata_10x_fp)

    # print out the saved adata
    print(adata_10x_subset, file=log_10x)
    print(adata_10x_subset.var, file=log_10x) 
    print(adata_10x_subset.obs, file=log_10x)       

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset
    level1_cts = adata_10x_subset.obs["cell_type_level_1"].dropna().unique()
    for ct in level1_cts:
        ct_data = adata_10x_subset[adata_10x_subset.obs["cell_type_level_1"] == ct, :].X
        out = ["level_1", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_10x)

    level2_cts = adata_10x_subset.obs["cell_type_level_2"].dropna().unique()
    for ct in level2_cts:
        ct_data = adata_10x_subset[adata_10x_subset.obs["cell_type_level_2"] == ct, :].X
        out = ["level_2", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_10x)

    level3_cts = adata_10x_subset.obs["cell_type_level_3"].dropna().unique()
    for ct in level3_cts:
        ct_data = adata_10x_subset[adata_10x_subset.obs["cell_type_level_3"] == ct, :].X
        out = ["level_3", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3_10x)

    table2_10x.close()
    table3_10x.close()
    log_10x.close()

main()


