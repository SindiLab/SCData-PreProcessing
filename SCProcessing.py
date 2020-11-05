import os
import numpy
import scanpy
import argparse
print("***SCProcessing Imported Successfully***")

# logging
import logging
from datetime import datetime
from Logging import LoggingHandler
from Preprocessing import Preprocess


logging.basicConfig(format='%(asctime)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S',
                    level=logging.INFO,
                    handlers=[LoggingHandler()])


# parse input arguments for preprocessing data
parser = argparse.ArgumentParser()
parser.add_argument('--minGenes', type=int, default=10, help='minimum number of genes, default=10')
parser.add_argument('--minCells', type=int, default=3, help='minimum number of cells, default=3')
parser.add_argument("--res", type=float, default=0, help="Clustering res, Default=0.15")
parser.add_argument("--randSeed", type=int, default=13, help="Random seed, Default=13")
parser.add_argument("--scaleFactor", type=int, default=20000, help="Random seed, Default=20000")
parser.add_argument("--rawDataPath", type=str, default=None, help="The path to raw data, this must be provided or will cause an error")
parser.add_argument("--savePath", type=int, default=None, help="The path to save the modified data and parameters, optional, Default=None")

global opt
opt = parser.parse_args()
print(opt)

if opt.savePath:
    os.makedirs(opt.savePath)

opt.rawDataPath =  "/home/jovyan/GitHub/DataSets/matrix.h5ad";
obj = Preprocess(opt.minGenes, opt.minCells, opt.res, opt.randSeed, opt.scaleFactor, opt.rawDataPath);
obj.Process();
