import random
import numpy as np
import pandas as pd
from natsort import natsorted
from collections import Counter, namedtuple

print("Coming to Train Split")

"""
To be added to a different file

#             'train_count': self.train_cells,
#             'valid_count': self.valid_cells,
#             'test_count': self.test_cells,
#             'train_cells_per_cluster': self.train_cells_per_cluster,
#             'valid_cells_per_cluster': self.valid_cells_per_cluster,
#             'test_cells_per_cluster': self.test_cells_per_cluster,
#             'clusters_no': self.clusters_no,
#             'clusters_ratios': self.clusters_ratios


"""



class TrainSplit():
    """
    
    A class for splitting the preprocessed data into train, validation and test splits 
    
    """
    
    def __init__(self, data, trainPercentage, validationPercentage, testPercentage, balancedSplit,randSeed):
        """
        
        Constructor of the TrainSplit class. It initializes class variables to be used in 
        splitting the data to train, test, validataion.
        
        """
        self.sc_raw = data
        self.train_cells = trainPercentage
        self.valid_cells = validationPercentage
        self.test_cells = testPercentage
        self.balanced_split = balancedSplit
        self.split_seed = randSeed
        self.cluster_res = clusterRes
        self.cells_count = self.sc_raw.shape[0]
        self.genes_count = self.sc_raw.shape[1]


    def Split(self):
        """
        Splits the data into training, validation and test sets.
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
            sc.tl.leiden(clustered, resolution = self.cluster_res, random_state = 1786)
#             sc.tl.louvain(clustered, resolution=)

            # add clusters to the raw data
#             self.sc_raw.obs['cluster'] = clustered.obs['louvain']
            self.sc_raw.obs['cluster'] = clustered.obs['leiden']


        # adding clusters' ratios
        cells_per_cluster = Counter(self.sc_raw.obs['cluster'])
        clust_ratios = dict()
        for key, value in cells_per_cluster.items():
            clust_ratios[key] = value / self.sc_raw.shape[0]

        self.clusters_ratios = clust_ratios
        self.clusters_no = len(cells_per_cluster)
        print("Clustering of the raw data is done to %d clusters."
              % self.clusters_no)

