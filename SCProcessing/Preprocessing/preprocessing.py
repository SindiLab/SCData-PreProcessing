import os 
import json
import logging
import warnings
logging.warn("Preprocessing: importing scanpy (may take a second)");
import scanpy as sc
logging.warn("Preprocessing: Importing Done");
import scipy.sparse as s_sparse

# Logging 
import logging
from datetime import datetime
# to hide pandas warning (or any other overhead libs)
hide_warnings = False;
if hide_warnings:
    import warnings
    warnings.filterwarnings("ignore")
logging.getLogger().setLevel(logging.INFO)

class Preprocess():
    """
    
    A class for preprocessing of single cell RNA-seq data
    
    """
    def __init__(self, minGene, minCells, res, randSeed, scaleFactor, rawDataPath):
        """
        
        Constructor of Preprocess class
        
        global attribute return:
            everything that has been passed on as an argument will be saved as class variables/attributes
        
        """
        self.min_genes = minGene;
        self.min_cells = minCells;
        self.seed = randSeed;
        self.cluster_res = res;
        self.scale = scaleFactor;
        
        if not rawDataPath:
            raise ValueError("You must provide the path raw dataset")
        else:
            self.raw_path = rawDataPath;
        
        
        
    def ReadData(self, inp_format_ = None):
        
        """
        
        method for reading raw data from a path and turning it into a dense matrix
        
        inputes:
            (optional) file format
        
        return:
            None
        
        class variable return:
            an anndata object from the read single cell data
            a format type for saving the file (later in Save method)
            a potential path for saving 
            
        modification:
            None
        
        """
        logging.info(f"Reading raw data from {self.raw_path}");
        
        # to see which file format has been provided to the method
        if inp_format_:
            self.inp_format = inp_format_
        else: 
            # split the string based on the '.' (because format comes right after)
            path = self.raw_path.split('.')
            # whatever comes after '.'
            self.inp_format = path.pop();
            # potential path for saving the modified data (in the same directory as the original data)
            self.save_path = path[0];
        
        # checking the different input formats 
        if self.inp_format == 'h5ad':
            adata = sc.read_h5ad(self.raw_path)
        
        elif self.inp_format == 'h5':
            adata = sc.read_10x_h5(self.raw_path)
        
        elif self.inp_format == 'mtx':
            adata = scanpy.read_mtx(self.raw_path)
        
        else:
            try:
                adata = scanpy.read(self.raw_path)
                
            except:
                raise ValueError(f"We have not gotten to {self.inp_format} yet, please explicitly specify this format for the method call")

        # keep track of unique gene names 
        adata.var_names_make_unique()
        
        if s_sparse.issparse(adata.X):
            print("Input data is already sparse");
            adata.X = adata.X.toarray()
        
        self.sc_raw = adata;
        
        logging.info(f"Read Data: {self.sc_raw}");
        
        
    def Filter(self):
        """
        Filtering Single cell data:
        
        inputes:
            None
        
        returns:
            None
        
        class variable return:
            new number of cells 
            new number of cells
            
        modification:
            it will directly modify the raw data that has been read
            
        """
        
        logging.info(f"BEFORE FILTERING: # of Cells: {self.sc_raw.shape[0]} -- # of Genes {self.sc_raw.shape[1]}");
        logging.info(f"Filtering cells less than {self.min_genes}")
        sc.pp.filter_cells(self.sc_raw, min_genes=self.min_genes, copy=False)
        
        logging.info(f"Filtering genes less than {self.min_cells}")
        sc.pp.filter_genes(self.sc_raw, min_cells=self.min_cells, copy=False)
        
        # record for later use
        self.cells_count = self.sc_raw.shape[0]
        self.genes_count = self.sc_raw.shape[1]
        
        logging.info(f"AFTER FILTERING: # of Cells: {self.cells_count} -- # of Genes {self.genes_count}")
        
        
        
    def Normalize(self):
        
        """
        Normalize library size given a factor 
        
        *** OTHER METHOD SHOULD BE IMPLEMENTED HERE IN THE FUTURE ***
        
        inputes:
            None
        
        returns:
            None
        
        class variable return:
            None
            
        modification:
            it will directly modify the raw data based on the scaled version
        
        """

        logging.info(f"Normalizing the data with a library size of {self.scale}")
        sc.pp.normalize_per_cell(self.sc_raw, counts_per_cell_after=self.scale)


    def Save(self, path = None):
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
        
        if path:
            self.save_path = path;
        else:
            str2remove = self.save_path.split("/").pop();
            self.save_path = self.save_path.replace(str2remove,'');
            
        save_path = self.save_path + str2remove + "_processed." + self.inp_format;
        logging.info(f"Saving data and parameters to folder {save_path}");
        print(f"SAVE PATH {save_path}")
        
        # write out the modified scanpy object
        self.sc_raw.write(save_path)
        
        # save all the parameters in a json file for reproducibility
        hparam = dict();
        
        hparam['preprocessed'] = \
        {
            'total_count': self.cells_count,
            'genes_no': self.genes_count,
            'split_seed': self.seed,
            'scale': self.scale
        }        
        
        with open(os.path.join(self.save_path, 'preprocessing_parameters.json'), 'w') as fp:
#             json.dump(hparam, fp, sort_keys=True, indent=4)
            json.dump(hparam, fp)
        
        
    def Process(self, save = True):
        """

        To automate the preprocessing procedure in one go!

        returns:
            None

        class variable return:
            method based

        modification:
            method based

        """
        
        self.ReadData();
        self.Filter();
        self.Normalize();
        if save:
            self.Save();
        else:
            logging.info("Skipping save as it was passed to this method")
        
        
        
    def GetAnnData(self):
        """
        
        Getter method for the modified preprocessed data
        
        returns:
            the annotated data
        
        class variable return:
            None

        modification:
            None
        
        """
        
        return self.sc_raw