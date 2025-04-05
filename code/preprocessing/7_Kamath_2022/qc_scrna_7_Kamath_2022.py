# This script processes the scRNAseq data from 7_Kamath2022_SN
# This script is run after running the script code\preprocessing\7_Kamath_2022\subset_for_normal.py
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2024-05-30

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import os
import single_cell_helper_functions_v3
from scipy.sparse import csc_matrix
import sys

# usage: python qc_scrna_7_Kamath_2022.py

def main():
# setting up #####
    base_dir = "Preprocessing_scRNA/data" #change here if necessary
    id = "426_Kamath2022_Human_2022_SN"

    # first concatenate
    opc = anndata.read("Preprocessing_scRNA/data/Kamathetal2022_SN/opc_normal.h5ad")
    endothelial = anndata.read("Preprocessing_scRNA/data/Kamathetal2022_SN/endothelial_normal.h5ad")
    microglia = anndata.read("Preprocessing_scRNA/data/Kamathetal2022_SN/microglia_normal.h5ad")
    astrocyte = anndata.read("Preprocessing_scRNA/data/Kamathetal2022_SN/astrocytes_normal.h5ad")
    da_neurons = anndata.read("Preprocessing_scRNA/data/Kamathetal2022_SN/DA_neurons_normal.h5ad")
    nurr_negative = anndata.read("Preprocessing_scRNA/data/Kamathetal2022_SN/nurr_negative_normal.h5ad")
    nurr_positive = anndata.read("Preprocessing_scRNA/data/Kamathetal2022_SN/nurr_positive_normal.h5ad")
    oligodendrocyte = anndata.read("Preprocessing_scRNA/data/Kamathetal2022_SN/oligodendrocytes_normal.h5ad")
    non_da_neurons = anndata.read("Preprocessing_scRNA/data/Kamathetal2022_SN/non_DA_neurons_normal.h5ad")
    
    adata = anndata.concat([opc, endothelial, microglia, astrocyte, da_neurons, nurr_negative, nurr_positive, oligodendrocyte, non_da_neurons])

    adata.var = opc.var

    gene_names_org = "symbol" #change here per data

    #### format output files ####
    log = open(os.path.join(base_dir, id, "log.txt"), "w")
    clean_adata = os.path.join(base_dir, id, id + ".h5ad")
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
    adata.obs.rename(columns={"author_cell_type": "cell_type_level_2"}, inplace = True)

    # modify the adata.var so that we have the row index of adata.var is symbol and ensemble column is named ensembl
    adata.var.rename(columns={"ensg_id": "ensg_id_orig"}, inplace = True) #the row name ensg_id exists both as row name and column name so I'm renaming the column here
    adata.var.reset_index(inplace=True) #since row index of adata.var is ensemble_id, we want to convert this into a column
    adata.var = adata.var.rename(columns={'feature_name': 'symbol'}) #update the column names to have consistent naming convention
    # convert to ensembl from symbol
    gene_names_fp = 'Preprocessing_scRNA/code/conversion_files/gene_names_human.txt' #path to gene_names_human.txt
    gene_names_add = "Ensembl gene ID"
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log)
    adata.var.set_index('symbol', inplace=True) # set the index to gene symbol
    # save to table 2
    ensembl_convertable = ["Ngenes with ensembl", str(adata.n_obs), str(ngenes_after_conversion)]
    print(",".join(ensembl_convertable), file=table2)

    # save
    adata.write_h5ad(filename=clean_adata)

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

    level2_cts = adata.obs["cell_type_level_2"].dropna().unique()
    for ct in level2_cts:
        ct_data = adata[adata.obs["cell_type_level_2"] == ct, :].X
        out = ["level_2", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)

    table2.close()
    table3.close()
    log.close()

main()