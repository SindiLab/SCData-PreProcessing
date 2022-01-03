import os
import json
import time
import random
import logging 
import numpy as np
import pandas as pd
from natsort import natsorted
print(" ")
logging.warn("TrainingSplit: importing scanpy (may take a second)");
import scanpy as sc
logging.warn("TrainingSplit: importing Done");
print(" ")

from collections import Counter, namedtuple
# to hide pandas warning (or any other overhead libs)
hide_warnings = False;
if hide_warnings:
    import warnings
    warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.INFO)


class TrainSplit():
    """
    
    A class for splitting the preprocessed data into train, validation and test splits 
    
    """
    
    def __init__(self, data, trainNumber, validationNumber, testNumber, balancedSplit:bool=True, randSeed:int=0, clusterRes=None, savePath = None, pre_process=False):
        
        """
        
        Constructor of the TrainSplit class. It initializes class variables to be used in 
        splitting the data to train, test, validataion.
        
        """
        self.sc_raw = data;
        self.train_cells = trainNumber;
        self.valid_cells = validationNumber;
        self.test_cells = testNumber;
        self.balanced_split = balancedSplit;
        self.split_seed = randSeed;
        self.cluster_res = clusterRes;
        self.cells_count = self.sc_raw.shape[0];
        self.genes_count = self.sc_raw.shape[1];
        self.save_path = savePath;
        self.pre_process = pre_process;

    def Split(self):
        
        """
        Class method for splitting the data into training, validation and test sets.
        If we want balanced splitting, this method will ensure the cluster ratios 
        are respected in each split.
        
        inputes:
            None 
        
        returns:
            None
        
        class variable return:
            None
            
        modification:
            None
        
        
        """
        print("==> Splitting:")
        start = time.time()
        random.seed(self.split_seed)
        np.random.seed(self.split_seed)

        if self.balanced_split:
            
            if hasattr(self, 'clusters_ratios'):
                logging.info("    -> Cluster Ratios exist")
                
                logging.info("    -> Starting a *balanced* split")
                print("    -> Starting a *balanced* split")

                valid_cells_per_cluster = {
                    key: int(value * self.valid_cells)
                    for key, value in self.clusters_ratios.items()}

                test_cells_per_cluster = {
                    key: int(value * self.test_cells)
                    for key, value in self.clusters_ratios.items()}

                dataset = np.repeat('train', self.sc_raw.shape[0])
                unique_groups = np.asarray(['valid', 'test', 'train'])
                self.sc_raw.obs['split'] = pd.Categorical(
                    values=dataset,
                    categories=natsorted(unique_groups))

                for key in valid_cells_per_cluster:

                    # all cells from clus idx
                    try:
                        indices = self.sc_raw.obs[
                            self.sc_raw.obs['cluster'] == str(key)].index
                    except:
                        logging.info("    -> Could not find the attribute 'cluster' in data")
                        logging.info("    -> Make sure this exists first before split")

                        
                    test_valid_indices = np.random.choice(
                        indices, valid_cells_per_cluster[key] +
                        test_cells_per_cluster[key], replace=False)

                    test_indices = test_valid_indices[0:test_cells_per_cluster[key]]
                    valid_indices = test_valid_indices[test_cells_per_cluster[key]:]

                    for i in test_indices:
                        self.sc_raw.obs.at[i, 'split'] =  'test'

                    for i in valid_indices:
                        self.sc_raw.obs.at[i, 'split'] = 'valid'

                self.valid_cells_per_cluster = valid_cells_per_cluster
                self.test_cells_per_cluster = test_cells_per_cluster
                
            else:
                print("GOT HERE")
                try: 
                    self.Cluster_ratios()
                except:
                    raise ValueError(" Cluster Ratios must exist... run object.Cluster() first")

        else:

            logging.info("    -> Starting a *non-balanced* split")
            print("    -> Starting a NON-balanced split")
            
            dataset = np.repeat('train', self.sc_raw.shape[0])
            unique_groups = np.asarray(['valid', 'test', 'train'])

            self.sc_raw.obs['split'] = pd.Categorical(
                values=dataset,
                categories=natsorted(unique_groups))

            # all cells from clus idx
            indices = self.sc_raw.obs.index

            test_valid_indices = np.random.choice(
                indices,
                self.test_cells + self.valid_cells,
                replace=False)

            test_indices = test_valid_indices[0:self.test_cells]
            valid_indices = test_valid_indices[self.test_cells:]

            for i in test_indices:
                self.sc_raw.obs.at[i, 'split'] =  'test'

            for i in valid_indices:
                self.sc_raw.obs.at[i, 'split'] = 'valid'

            self.valid_cells_per_cluster = Counter(
                self.sc_raw[valid_indices].obs['cluster'])

            self.test_cells_per_cluster = Counter(
                self.sc_raw[test_indices].obs['cluster'])

        train_indices = self.sc_raw[
            self.sc_raw.obs['split'] == 'train'].obs.index

        self.train_cells_per_cluster = dict(
            Counter(self.sc_raw[train_indices].obs['cluster']))

        self.train_cells = self.cells_count - self.test_cells - self.valid_cells
        
        print("-><- Splitting done")
        print(f"Splitting took {time.time() - start} seconds")
        
    def Cluster(self, clusterRes=None):
        
        """
        Method that applies a Louvain clustering of the data, following
        the Zheng recipe. Computes and stores the cluster ratios.
        
        inputes:
            None 
        
        returns:
            None
        
        class variable return:
            None
            
        modification:
            None
        
        """
        print("==> Starting to cluster")
        start = time.time()
        if self.cluster_res is None:
            if "cluster" in self.sc_raw.obs_keys():
                print("    -> clustering is already done, no clustering will be applied.")
                print("-><- Returning")
                # returning
                return 0
        
        if clusterRes is None:
            if self.cluster_res is None:
                raise ValueError('please provide a resolution (i.e. set clusterRes) and run Cluster() again :) ')
        else:
            print(f"    -> Setting clustering resolution to {clusterRes}")
            self.cluster_res = clusterRes
       
        clustered = self.sc_raw.copy()

        # pre-processing
        if self.pre_process:
            print("    -> Pre-processing (as set True by user):")
            sc.pp.recipe_zheng17(clustered)
            print("    -> Pre-processing done.")
        
        if "X_pca" in self.sc_raw.obsm_keys():
            print("    -> Found existing PCA reduction within object. Skipping PCA")
        else:
            print("    -> Running PCA:")
            sc.tl.pca(clustered, n_comps=50)
            self.sc_raw.obsm['X_pca'] = clustered.obsm['X_pca']
            self.sc_raw.varm['PCs'] = clustered.varm['PCs']
            print("    -> PCA done.")
            
        # clustering
        print("    -> Clustering:")
        sc.pp.neighbors(clustered, n_pcs=50)
        sc.tl.leiden(clustered, resolution = self.cluster_res, random_state = self.split_seed)
        """
        ###### depreciated
        sc.tl.louvain(clustered, resolution=)
        # add clusters to the raw data
        self.sc_raw.obs['cluster'] = clustered.obs['louvain']
        """
        self.sc_raw.obs['cluster'] = clustered.obs['leiden']
        # get cluster ratios for balances split
        self.Cluster_ratios();
        print("-><- Done. Clustering of the raw data is done to %d clusters." % self.clusters_no)
        print(f" Clustering took {time.time() - start} seconds")

        
    def Cluster_ratios(self):
        print("==> Saving cluster ratios:")
        # adding clusters' ratios
        try:
            cells_per_cluster = Counter(self.sc_raw.obs['cluster'])
        except:
            self.sc_raw.obs['cluster'] = self.sc_raw.obs['louvain']
            cells_per_cluster = Counter(self.sc_raw.obs['cluster'])
            
        clust_ratios = dict()
        for key, value in cells_per_cluster.items():
            clust_ratios[key] = value / self.sc_raw.shape[0]

        self.clusters_ratios = clust_ratios
        self.clusters_no = len(cells_per_cluster)
        print(f"    -> Number of clusters: {self.clusters_no}")
        print("-><- Saved cluster ratios to object attributes")
        
        
    def Save(self, save_path = None):
        
        """
        Save the modified data and the associated parameters
        
        inputes:
            (optional) path: path to the directory we want to save the data+params in 
        
        returns:
            None
        
        class variable return:
            None
            
        modification:
            None
        
        """
        print("==> Saving processed data:")
        if save_path != None:
            self.save_path = save_path;
            
        if not self.save_path:
            # let's try to make a new directory
            try:
                os.makedirs("./TrainTestSplitData")
            except:
                pass
            
            self.save_path = "./TrainTestSplitData/"
            
        else:
            print(f"    -> trying to save in: {self.save_path}")
            try:
                os.makedirs(f"{self.save_path}")
            except:
                pass
            
        save_path = self.save_path + "TrainTestSplit" + ".h5ad";
        
        print(f"    -> Saving data and parameters to folder {save_path}");
        
        # write out the modified scanpy object
        self.sc_raw.write(save_path)
        
        # save all the parameters in a json file for reproducibility
        hparam = dict();

        hparam['preprocessed'] = \
        {
            'clustering_resolution': self.cluster_res,
            'train_count': self.train_cells,
            'valid_count': self.valid_cells,
            'test_count': self.test_cells,
            'train_cells_per_cluster': self.train_cells_per_cluster,
            'valid_cells_per_cluster': self.valid_cells_per_cluster,
            'test_cells_per_cluster': self.test_cells_per_cluster,
            'clusters_no': self.clusters_no,
            'clusters_ratios': self.clusters_ratios
        }        
        
        with open(os.path.join(self.save_path, 'TrainTestSplit_parameters.json'), 'w') as fp:
#             json.dump(hparam, fp, sort_keys=True, indent=4)
            json.dump(hparam, fp)
    
        print("-><- Saving done.")
        
    def GetData(self):
        
        return self.sc_raw;