
import argparse
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn import model_selection
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import cross_val_score
from sklearn.svm import SVC
from scipy import sparse
import os
import json
import tszip
import matplotlib.pyplot as plt
import seaborn as sns
import pickle

def get_1kg_genotypes(snp_df, tree_seq_dir):
    chroms = snp_df.chrom.unique()
    snp_df['merge'] = [f"{s.chrom}_{s.pos}" for i,s in snp_df.iterrows()]
    
    chr_dfs = []
    genotype_matrices = []
    
    for chrom in chroms:
        ts_path = os.path.join(tree_seq_dir, f"1kg_chr{chrom}.trees.tsz")
        if os.path.isfile(ts_path): 
            t0 = time()
            chr_ts = tszip.decompress(ts_path)
            muts = pd.DataFrame([{'merge': f"{chrom}_{int(m.position)}",
                                  'ts_id': m.site,
                                  'alt':m.derived_state} for m in chr_ts.mutations()])
            shared_muts = snp_df.merge(muts, how='inner', on=['merge']).drop_duplicates(subset='ts_id')
            unshared_set = list(set(muts.ts_id) - set(shared_muts.ts_id))
            chr_ts_shared = chr_ts.delete_sites(unshared_set)
            g_1kg_haploid = chr_ts_shared.genotype_matrix().T
            g_1kg = g_1kg_haploid[0::2] + g_1kg_haploid[1::2]
            
            chr_dfs.append(shared_muts)
            genotype_matrices.append(sparse.csr_matrix(g_1kg))
            t1 = time() - t0
            print(f"Got genotypes for chromosome {chrom} with with shape {g_1kg.shape} in {t1:.2f} seconds")
            
        else:
            print(f"No file {ts_path.split('/')[-1]} for chromosome {chrom}")

    genotype_matrix_1kg = sparse.hstack(genotype_matrices)
    variant_df_1kg = pd.concat(chr_dfs)
    variant_df_1kg.drop(columns=['merge', 'ts_id'])
    return genotype_matrix_1kg, variant_df_1kg

parse_metadata = lambda metadata: json.loads(metadata.decode())

def get_1kg_population_df(tree_seq_dir):
    """ Read in population and subpopulation for individuals """
    chrom = 22
    ts_path = os.path.join(tree_seq_dir, f"1kg_chr{chrom}.trees.tsz")
    chr_ts = tszip.decompress(ts_path)

    individual_ids = [parse_metadata(individual.metadata)['individual_id']
                      for individual in chr_ts.individuals()]
    populations = pd.DataFrame([parse_metadata(population.metadata)
                                for population in chr_ts.populations()])
    num_individuals = chr_ts.num_individuals
    individual_population_ids = [chr_ts.node(2*i).population for i in range(num_individuals)]

    individual_populations = list(populations.name[individual_population_ids])
    individual_super_populations = list(populations.super_population[individual_population_ids])

    pop_df = pd.DataFrame({'id':individual_population_ids,
                         'pop':individual_populations,
                         'superpop':individual_super_populations})

    return pop_df


def plot_confusion(Y_pred, Y_true, label_encoder, ax=None):
    labels = label_encoder.classes_
    Y_pred_l = label_encoder.inverse_transform(Y_pred)
    Y_true_l = label_encoder.inverse_transform(Y_true)
    C = confusion_matrix(Y_pred_l, Y_true_l, labels=labels)

    if ax:
        sns.heatmap(C, annot=True, fmt='d', ax=ax)
    else:
        ax = plt.subplot()
        sns.heatmap(C, annot=True, fmt='d')

    ax.xaxis.set_ticklabels(labels)
    ax.yaxis.set_ticklabels(labels)
    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")
    return


def train_svm(G, pop_df, aim_df):
    """ This trains the SVM on 1KG genotypes 
    and evaluates the performance with 10-fold cross validation """

    # Turn populations into numbers
    pops = pop_df['pop'].unique()
    superpops = pop_df['superpop'].unique()

    le_pop = preprocessing.LabelEncoder().fit(pops)
    le_superpop = preprocessing.LabelEncoder().fit(superpops)

    Y_super = le_superpop.transform(pop_df['superpop'])
    Y_pop = le_pop.transform(pop_df['pop'])
    
    # Initial test: 90/10 train test split
    # Get random train/test split
    G = G.todense()
    train_ind_rand, test_ind_rand = next( 
        model_selection.ShuffleSplit(
            n_splits=1,
            train_size=0.9, 
            test_size=0.1).split(G, Y_super))

    G_train, Y_train = G[train_ind_rand], Y_super[train_ind_rand]
    G_test, Y_test = G[test_ind_rand], Y_super[test_ind_rand]

    clf = SVC().fit(G_train, Y_train)
    Y_test_pred = clf.predict(G_test)

    plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plot_confusion(Y_test_pred, Y_test, le_superpop, ax)
    ax.collections[0].colorbar.set_label("Num. of individuals")
    plt.title("Test performance on 1KG")
    plt.tight_layout()
    plt.savefig("SVC_perf_on_1KG.png", dpi=400)
    plt.clf()

    # Subsequent test: 10-fold cross validation
    scores = cross_val_score(SVC(), G, Y_super, cv=10)
    plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    plt.title("10-fold cross-validated accuracy on 1KG")
    sns.distplot(scores, kde=False)
    plt.xlabel("Accuracy")
    plt.ylabel("Number of folds")
    plt.savefig("SVC_10fold_cv_accuracy.png", dpi=400)

    # Now, train on full dataset and return classifier
    clf_full = SVC().fit(G, Y_super)
    clf_full.snp_names = aim_df.rsid.values
    return clf_full, le_superpop


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

    aim_df = pd.read_csv(args.input, sep=" ")

    precomp_genotype_file = os.path.join(args.data, "genotype_matrix.npz")
    precomp_aim_file = os.path.join(args.data, "aim_variants.pkl")
    # This is the slow step
    if (os.path.isfile(precomp_genotype_file) and os.path.isfile(precomp_aim_file)):
       print("True")
       aim_G = sparse.load_npz(precomp_genotype_file)
       aim_df = pd.read_pickle(precomp_aim_file)
    else:
        print("False")
        aim_G, aim_df = get_1kg_genotypes(args.data)

    pop_df =  get_1kg_population_df(args.data)
    clf, le = train_svm(aim_G, pop_df, aim_df)
    
    trained_clf_file = os.path.join(args.data, "ancestry_svc_clf.pkl")
    with open(trained_clf_file, "wb") as f:
        f.write(pickle.dumps(clf))

    trained_le_file = os.path.join(args.data, "ancestry_label_encoder.pkl")
    with open(trained_le_file, "wb") as f:
        f.write(pickle.dumps(le))

    return 

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", 
        help="Path to text file containing the AIMs used for training",
        default="data/kidd_et_al_aims.txt")
    parser.add_argument("-d", "--data", 
        help="Path to directory containing the 1000 genomes tree seqences",
        default="data")

    main(parser.parse_args())
