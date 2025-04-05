# This script subset for normal cells only
# The rationale is that when I downloaded the h5ad file for the Kamath dataset, it's divided into each cell type. Therefore, the idea is that I will subset to keep the normal cells only and then will concatenate and then will filter on the full file
# Author: Tanya Phung (t.n.phung@vu.nl)
# Date: 2024-05-30

import scanpy as sc
import pandas as pd
import numpy as np
import anndata
import os
from scipy.sparse import csc_matrix
import sys
# usage: python subset_for_normal.py h5ad_path out_h5ad_path
def main():
# setting up #####
    h5ad_path = sys.argv[1]
    normal_h5ad_path = sys.argv[2]

    #### read in adata #####
    adata = anndata.read(h5ad_path)

    metadata = adata.obs
    metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
    metadata = metadata.rename(columns = {'index': 'cell_id'})

    # separate for normal only
    print("Subsetting the adata for normal")

    cell_id = metadata[metadata["disease"] == "normal"]["cell_id"]
    adata.obs.set_index('index', inplace=True)
    adata_subset = adata[cell_id]

    print("Subsetting adata for group is finished.")

    # save
    adata_subset.write_h5ad(filename=normal_h5ad_path)

main()