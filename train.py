import argparse
import pandas as pd

def main(args):
    """
    Trains a model inferring genetic ancestry from a pre-specified
    set of common ancestry informative marker (AIMs). 

    Performs the following steps:
    1. Creates the training data from 1000 genomes sequences.
    This involves extracting the genotypes for each individual at each AIM
    along with their super-population information. Because we are working with
    and extremely compressed format of the 1000 genomes sequences, this step takes 
    ~20-30 minutes on a laptop. 

    2. Trains a support vector machine to predict super-population label from genotype.
    The model is saved and diagnostic outputs are produced.
    """
    return 

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", 
        help="Path to text file containing the AIMs used for training",
        default="data/kidd_et_al_aims.txt")
    parser.add_argument("-t", "--tree-sequences", 
        help="Path to directory containinf the 1000 genomes tree seqences",
        default="data")

    main(parser.parse_args())
