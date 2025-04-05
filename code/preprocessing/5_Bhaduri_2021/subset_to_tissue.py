import scanpy as sc
import pandas as pd
import numpy as np
import anndata
from scipy.sparse import csc_matrix
pd.set_option('display.max_columns', None)
pd.set_option('display.max_row', None)

# subset to 6 tissues to reduce size
print("Subsetting the adata for each tissue")

tissue_dict = {"PrefrontalCortex": "prefrontal cortex",
               "PrimaryMotorCortex": "primary motor cortex",
               "ParietalCortex": "parietal cortex",
               "PrimarySomatosensoryCortex": "primary somatosensory cortex",
               "TemporalCortex": "temporal cortex",
               "PrimaryVisualCortex": "primary visual cortex"}
for key in tissue_dict:
    adata = anndata.read("d0bbb292-7b5d-4f49-b3d8-046bc32a444c.h5ad")
    # make the X layer to have the raw count layer
    print("Make the X layer to have the raw count layer")
    adata.X = adata.raw.X.copy()
    print("Finished with making the X layer to have the raw count layer")
    metadata = adata.obs
    metadata.reset_index(inplace=True) #because the cell id is set as the row index, we convert the row index into a column
    metadata = metadata.rename(columns = {'index': 'cell_id'})
    cell_id = metadata[(metadata["tissue"] == tissue_dict[key])]["cell_id"]
    adata.obs.set_index('index', inplace=True)
    adata_subset = adata[cell_id]
    out = key + ".h5ad"
    adata_subset.write_h5ad(filename=out)
    print("Finish for ", key)