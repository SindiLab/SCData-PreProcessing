import os
import json
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
    
    def __init__(self, data, trainPercentage, validationPercentage, testPercentage, balancedSplit, randSeed, clusterRes=None, savePath = None):
        
        """
        
        Constructor of the TrainSplit class. It initializes class variables to be used in 
        splitting the data to train, test, validataion.
        
        """
        self.sc_raw = data;
        self.train_cells = trainPercentage;
        self.valid_cells = validationPercentage;
        self.test_cells = testPercentage;
        self.balanced_split = balancedSplit;
        self.split_seed = randSeed;
        self.cluster_res = clusterRes;
        self.cells_count = self.sc_raw.shape[0];
        self.genes_count = self.sc_raw.shape[1];
        self.save_path = savePath;

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
        
        random.seed(self.split_seed)
        np.random.seed(self.split_seed)

        if self.balanced_split:
            logging.info("Starting a *balanced* split")
            print("Starting a *balanced* split")

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
                indices = self.sc_raw.obs[
                    self.sc_raw.obs['cluster'] == str(key)].index

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

            logging.info("Starting a *non-balanced* split")
            print("Starting a NON-balanced split")
            
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
        
        logging.info("Splitting done")
        
        
        
    def Cluster(self):
        
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
        logging.info("Starting to cluster")
        
        if self.cluster_res is None:
            if "cluster" in self.sc_raw.obs_keys():
                print("clustering is already done,"
                      " no clustering will be applied")
            else:
                raise ValueError(' No clustering is applied, '
                                 'please apply clustering')
        else:
            clustered = self.sc_raw.copy()

            # pre-processing
            sc.pp.recipe_zheng17(clustered)
            sc.tl.pca(clustered, n_comps=50)

            # clustering
            sc.pp.neighbors(clustered, n_pcs=50)
            sc.tl.leiden(clustered, resolution = self.cluster_res, random_state = self.split_seed)
            """
            ###### depreciated
            sc.tl.louvain(clustered, resolution=)
            # add clusters to the raw data
            self.sc_raw.obs['cluster'] = clustered.obs['louvain']
            """
            self.sc_raw.obs['cluster'] = clustered.obs['leiden']


        # adding clusters' ratios
        cells_per_cluster = Counter(self.sc_raw.obs['cluster'])
        clust_ratios = dict()
        for key, value in cells_per_cluster.items():
            clust_ratios[key] = value / self.sc_raw.shape[0]

        self.clusters_ratios = clust_ratios
        self.clusters_no = len(cells_per_cluster)
        
        logging.info("Clustering of the raw data is done to %d clusters." % self.clusters_no)


        
    def Save(self):
        
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
        
        if not self.save_path:
            # let's try to make a new directory
            try:
                os.makedirs("./TrainSplitData")
            except:
                pass
            
            self.save_path = "./TrainSplitData/"
            
        save_path = self.save_path + "TrainSplit" + ".h5ad";
        
        logging.info(f"Saving data and parameters to folder {save_path}");
        
        # write out the modified scanpy object
        self.sc_raw.write(save_path)
        
        # save all the parameters in a json file for reproducibility
        hparam = dict();

        hparam['preprocessed'] = \
        {
            'train_count': self.train_cells,
            'valid_count': self.valid_cells,
            'test_count': self.test_cells,
            'train_cells_per_cluster': self.train_cells_per_cluster,
            'valid_cells_per_cluster': self.valid_cells_per_cluster,
            'test_cells_per_cluster': self.test_cells_per_cluster,
            'clusters_no': self.clusters_no,
            'clusters_ratios': self.clusters_ratios
        }        
        
        with open(os.path.join(self.save_path, 'TrainSplit_parameters.json'), 'w') as fp:
#             json.dump(hparam, fp, sort_keys=True, indent=4)
            json.dump(hparam, fp)