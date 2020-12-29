import os
import argparse
print("*** SCProcessing Imported Successfully ***")
# logging 
import logging
from datetime import datetime
# Package classes
from SCProcessing import Preprocess, TrainSplit, LoggingHandler, read_and_serialize

import scanpy as sc

# if you don't want any logging output from the class methods
verbose = True;
if verbose:
    logging.basicConfig(format='%(asctime)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO,
                    handlers=[LoggingHandler()])


# parse input arguments for preprocessing data 
parser = argparse.ArgumentParser()
parser.add_argument('--RawDataDirectory', type=str, default=None, help='directory for the raw data, default=None')
parser.add_argument('--minGenes', type=int, default=10, help='minimum number of genes, default=10')
parser.add_argument('--minCells', type=int, default=3, help='minimum number of cells, default=3')
parser.add_argument("--scaleFactor", type=int, default=20000, help="Random seed, Default=20000")
parser.add_argument("--rawDataPath", type=str, default=None, help="The path to raw data, this must be provided or will cause an error")
parser.add_argument("--savePath", type=int, default=None, help="The path to save the modified data and parameters, optional, Default=None")

# parse input arguments for splitting the data into train, validation and test
# parser.add_argument("--trainPercentage", type=int, default=90, help="Percentage of training data, Default=80")
# parser.add_argument("--validationPercentage", type=int, default=10, help="Percentage of validation data, Default=20")

parser.add_argument("--validationCells", type=int, default=2000, help="Number of validation cells, Default=2000 for 68K")
parser.add_argument("--testPercentage", type=int, default=0, help="Percentage of testing data, Default=0")
parser.add_argument("--balancedSplit", type=bool, default=True, help="Whether we want a balance split or not, Default=True")
parser.add_argument("--randSeed", type=int, default=0, help="Random seed, Default=0")
parser.add_argument("--res", type=float, default=0.15, help="Clustering res, Default=0.15")




global opt
opt = parser.parse_args()
print(opt)

if not opt.RawDataDirectory:
    opt.rawDataPath =  "/home/jovyan/68K_PBMC_Raw/raw_68kPBMCs.h5ad";

if opt.savePath:
    os.makedirs(opt.savePath)
else:
    opt.savePath = "/home/jovyan/68K_PBMC_Processed/68kPBMCs_processedLikeSCGAN.h5ad";

# if we have to prepare the data from raw SC data

obj = Preprocess(opt.minGenes, opt.minCells, opt.res, opt.randSeed, opt.scaleFactor, opt.rawDataPath);
obj.Process();

# get the modified data so that we can do training splits and clustering 
proc_data = obj.GetAnnData();

train_split = TrainSplit(proc_data, opt.trainPercentage, opt.validationPercentage, opt.testPercentage, opt.balancedSplit, opt.randSeed, opt.res, opt.savePath);

print(proc_data);
train_split.Cluster()
train_split.Split()
train_split.Save()
logging.info("Train, Valid and Test split is complete")
split_data = train_split.GetData()


## if the data with the splits is already available ##
# split_data = sc.read("/home/jovyan/GitHub/DataSets/TrainSplit.h5ad")
print(split_data);


# print("Calling TF Creator")
# read_and_serialize(split_data,opt.savePath) 
print("DONE")