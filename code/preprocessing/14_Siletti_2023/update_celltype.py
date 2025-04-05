# This script updates `cell_type` to `cell_type_level_1` and `supercluster_term` to `cell_type_level_2`
# Author: Tanya Phung (t.n.phung@vu.nl)

import scanpy as sc
import pandas as pd
import anndata
import os
import argparse

def main(args):
    basedir = "Preprocessing_scRNA/data" #NOTE: change here if necessary
    id = args.id

    #rename
    original_fp = os.path.join(basedir, id, id + ".h5ad")
    temp_fp = os.path.join(basedir, id, id + "_tmp.h5ad")

    os.rename(original_fp, temp_fp)
    adata = anndata.read(temp_fp)
    # rename cell type columns

    adata.obs.rename(columns={"cell_type": "cell_type_level_1"}, inplace = True)
    adata.obs.rename(columns={"supercluster_term": "cell_type_level_2"}, inplace = True)
    
    # save
    adata.write_h5ad(filename=original_fp)

    # remove temp file
    os.remove(temp_fp)



def parse_args():
    parser = argparse.ArgumentParser(description='Update cell type columns in adata.obs to be consistent')
    parser.add_argument('--id', required=True, help='Insert ID. For example: 60_Siletti_Hippocampus.HiH.HiT.Sub_Human_2022')
    return parser.parse_args()


main(parse_args())


