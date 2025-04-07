# This script processes the scRNAseq data from Phan 2024
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2024-02-13

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import os
import single_cell_helper_functions_v3
from scipy.sparse import csc_matrix
import sys

# usage: python qc_scrna_19_Phan_2024.py 513_Pan_Human_2024_PrefrontalCortex PrefrontalCortex

def main():
# setting up #####
    id = sys.argv[1]
    tissue = sys.argv[2]
    
    tissue_dict = {
        "PrefrontalCortex": "prefrontal cortex",
        "WhiteMatterOccipitalLobe": "white matter of occipital lobe",
        "WhiteMatterTemporalLobe": "white matter of temporal lobe",
        "WhiteMatterParietalLobe": "white matter of parietal lobe",
        "WhiteMatterFrontalLobe": "white matter of frontal lobe"
    }

    base_dir = "Preprocessing_scRNA/data" #change here if necessary
    h5ad_path = "982952fb-dee0-4498-a072-3dc498b06aae.h5ad" #change here per data
    gene_names_org = "symbol" #change here per data
    
    # setting up output
    log = open(os.path.join(base_dir, id, "log.txt"), "w")
    clean_adata_fp = os.path.join(base_dir, id, id + ".h5ad")

    # format table 2
    table2 = open(os.path.join(base_dir, id, "table2.csv"), "w")
    header = ["Description", "Number of cells", "Number of genes"]
    print(",".join(header), file=table2)

    # format table 3
    table3 = open(os.path.join(base_dir, id, "table3.csv"), "w")

    # read in adata #####
    adata = anndata.read(h5ad_path)
    
    # subset for tissue and dev age
    metadata = adata.obs
    metadata.reset_index(inplace=True)
    cell_id = metadata[(metadata["tissue"] == tissue_dict[tissue]) & (metadata["disease"] == "normal")]["index"] 
    adata.obs.set_index('index', inplace=True)
    adata_subset = adata[cell_id]
    print("Subsetting adata for group is finished.", file=log)
    print("Viewing the adata observations.", file=log)
    print(adata_subset.obs, file=log)

    print("Viewing the adata variables.", file=log)
    print(adata_subset.var, file=log)
    
    print("Make the X layer to have the raw count layer", file=log)
    adata_subset.X = adata_subset.raw.X.copy()
    print("Finished with making the X layer to have the raw count layer", file=log)
    print("Viewing the adata matrix - are these integer counts?",file=log)
    print(adata_subset.X.A[1:25,1:25],file=log)

    adata_nrow = adata_subset.shape[0]
    adata_ncol = adata_subset.shape[1]
    print(f"After subsetting, adata has {adata_nrow} cells and {adata_ncol} genes.", file=log)
    original = ["Original", str(adata_subset.n_obs), str(adata_subset.n_vars)]
    print(",".join(original), file=table2)

    # mitochondrial genes
    adata_subset.var['mt'] = adata_subset.var["feature_name"].str.startswith("MT-")
    # count the number of genes that are mitochondrial genes
    n_mt_genes = np.count_nonzero(adata_subset.var['mt'])
    print(f"adata has {n_mt_genes} mitochondrial genes.", file=log)

    # calculate qc metrics
    sc.pp.calculate_qc_metrics(adata_subset, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # first filtering: accept cells with at least 200 detected genes and genes need to be expressed in at least 3 cells
    sc.pp.filter_cells(adata_subset, min_genes=200)
    sc.pp.filter_genes(adata_subset, min_cells=3)

    print(f"After first filtering (see documentation for definition), adata has {adata_subset.n_obs} cells and {adata_subset.n_vars} genes.", file=log)
    # save to table 2
    first_filter = ["Filter#1", str(adata_subset.n_obs), str(adata_subset.n_vars)]
    print(",".join(first_filter), file=table2)

    # second filtering based on mt percentage
    adata_subset = adata_subset[adata_subset.obs['pct_counts_mt'] < 10, :]
    print(f"After second filtering (see documentation for definition), adata has {adata_subset.n_obs} cells and {adata_subset.n_vars} genes.", file=log)
    # save to table 2
    second_filter = ["Filter#2", str(adata_subset.n_obs), str(adata_subset.n_vars)]
    print(",".join(second_filter), file=table2)

    # rename cell type columns
    adata_subset.obs.rename(columns={"cell_type": "cell_type_level_1"}, inplace = True)

    # modify the adata.var so that we have the row index of adata.var is symbol and ensemble column is named ensembl
    adata_subset.var.reset_index(inplace=True) #since row index of adata.var is ensemble_id, we want to convert this into a column
    adata_subset.var = adata_subset.var.rename(columns={'feature_name': 'symbol'}) #update the column names to have consistent naming convention
    # convert to ensembl from symbol
    gene_names_fp = 'gene_names_human.txt' #path to gene_names_human.txt
    gene_names_add = "Ensembl gene ID"
    adata_gene_converted = single_cell_helper_functions_v3.add_gene_names_human(adata_subset, gene_names_fp, gene_names_org, gene_names_add)
    ngenes_after_conversion = adata_gene_converted.var.dropna(subset=["ensembl"]).shape[0]
    print(f"After converting to ensembl, adata has {ngenes_after_conversion} genes.", file=log)
    adata_subset.var.set_index('symbol', inplace=True) # set the index to gene symbol
    # save to table 2
    ensembl_convertable = ["Ngenes with ensembl", str(adata_subset.n_obs), str(ngenes_after_conversion)]
    print(",".join(ensembl_convertable), file=table2)

    # save
    adata_subset.write_h5ad(filename=clean_adata_fp)

    # print out the saved adata
    print(adata_subset, file=log)
    print(adata_subset.var, file=log) 
    print(adata_subset.obs, file=log)       

    # tabulate the number of cells per cell type
    # define the available cell types in the dataset
    level1_cts = adata_subset.obs["cell_type_level_1"].dropna().unique()
    for ct in level1_cts:
        ct_data = adata_subset[adata_subset.obs["cell_type_level_1"] == ct, :].X
        out = ["level_1", ct, str(ct_data.shape[0])]
        print("|".join(out), file=table3)

    table2.close()
    table3.close()
    log.close()

main()