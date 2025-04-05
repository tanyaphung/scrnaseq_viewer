# This script preproprocesses and outputs qc metrics for all of the dissection for Silleti Adult Human Brain 2022
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2023-05-13
# Note that I have removed all paths and just keep the basename. For rerunning one would need to make sure that paths, etc... are set correctly. Search for keyword NOTE to see all the places one may need to change

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from matplotlib import pyplot as plt
import os
import single_cell_helper_functions_v3 #NOTE: need to make sure that the path is set correctly to load this helper func
import argparse

def main(args):
    id = args.id
    ct_colname = "supercluster_term"
    base_dir = "Preprocessing_scRNA/data" #NOTE: change here if necessary

    h5ad_path = os.path.join(base_dir, id, "local.h5ad") #NOTE: change here if necessary

    gene_names_org = "symbol" #NOTE: change here if necessary

    log = open(os.path.join(base_dir, id, "log.txt"), "w")
    plot1 = os.path.join(base_dir, id, "plot1.jpeg")
    data_for_plot2_fn = os.path.join(base_dir, id, id + "_obs_for_plot2.csv")

    clean_adata_fp = os.path.join(base_dir, id, id + ".h5ad")

    # format table 2
    table2 = open(os.path.join(base_dir, id, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2)

    # format table 3
    table3 = open(os.path.join(base_dir, id, "table3.csv"), "w")

    # read in adata
    adata = anndata.read(h5ad_path)

    print("Viewing the adata observations.", file=log)
    print(adata.obs, file=log)

    print("Viewing the adata variables.", file=log)
    print(adata.var, file=log)

    adata_nrow = adata.shape[0]
    adata_ncol = adata.shape[1]
    print(f"adata has {adata_nrow} cells and {adata_ncol} genes.", file=log)
    original = ["Original", str(adata.n_obs), str(adata.n_vars)]
    print(",".join(original), file=table2)

    # mitochondrial genes
    adata.var['mt'] = adata.var["Gene"].str.startswith("MT-")
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

    # modify the adata.var so that we have the row index of adata.var is symbol and ensemble column is named ensembl
    adata.var.reset_index(inplace=True) #since row index of adata.var is ensemble_ids, we want to convert this into a column
    adata.var = adata.var.rename(columns={'Gene': 'symbol'}) #update the column names to have consistent naming convention
    # convert to ensembl from symbol
    gene_names_fp = 'gene_names_human.txt' #path to gene_names_human.txt #NOTE: make sure path is correct. One can find this file at `resources/gene_names_human.txt` on this repor
    gene_names_add = "Ensembl gene ID"
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log)
    adata.var.set_index('symbol', inplace=True) # set the index to gene symbol

    # save
    adata.write_h5ad(filename=clean_adata_fp)

    # plot highest expressed genes
    with plt.rc_context():  # Use this to set figure params like size and dpi
        sc.pl.highest_expr_genes(adata, n_top=20, show=False)
        plt.savefig(plot1, bbox_inches="tight")

    # save data to plot violin in R
    data_for_plot2 = adata.obs.reset_index()
    data_for_plot2.to_csv(data_for_plot2_fn, index=False)

    # print out the saved adata
    print(adata, file=log)
    print(adata.var, file=log) 
    print(adata.obs, file=log)       

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset
    cts = adata.obs[ct_colname].dropna().unique()
    for ct in cts:
        ct_data = adata[adata.obs[ct_colname] == ct, :].X
        out = [ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)

    table2.close()
    table3.close()
    log.close()


def parse_args():
    parser = argparse.ArgumentParser(description='Preprocess scRNAseq data from Siletti et al. 2022 for adult human brain.')
    parser.add_argument('--id', required=True, help='Insert ID. For example: 7_Siletti_CerebralCortex.PrCG.M1C_Human_2022')
    return parser.parse_args()


main(parse_args())


