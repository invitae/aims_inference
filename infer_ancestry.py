"""
Predicts ancestry from a fixed set of SNPs

Requires a very specifically-formatted input file for genotypes SNPs:
    1. Comma-separated values, with the first line as a header
    2. First column is named 'sample' 
    3. Sample values have the format SAMPLENAME_ETHNICITY
    4. Subsequent columns are rsids of SNPs (ex. rs3814134)
    5. Genotypes are labeled like: A/C

"""

import numpy as np
import pickle
import pandas as pd
import argparse
from sklearn.metrics import confusion_matrix
from sklearn import preprocessing
import seaborn as sns
import matplotlib.pyplot as plt

def encode_test_data(input_file, aim_dataframe_file, parse_ethnicity=True):
    vdf = pd.read_pickle(aim_dataframe_file)
    vdf.reset_index(inplace=True, drop=True)

    test_df = pd.read_csv(input_file)
    test_df.reset_index(inplace=True)
    if parse_ethnicity:
        test_df['ethnicity'] = \
            test_df['sample'].str.split("_").apply(lambda x: x[1]).astype("category")
    test_aims = [col for col in test_df.columns if col.startswith("rs")]
    
    G_test = []
    for snp in vdf.rsid.values:
        snp_col = test_df[snp]
        allele = vdf[vdf.rsid==snp].alt.values[0]
        G_test.append(np.array(test_df[snp].str.count(allele)))
    G_test = np.vstack(G_test).T
    return G_test, test_df

def load_classifier(clf_file, le_file):
    with open(clf_file, "rb") as f:
        clf = pickle.load(f)
    
    with open(le_file, "rb") as f:
        le = pickle.load(f)
    return clf, le


def test_confusion_matrix(sample_df, anc_le, output_filename):
    ethnicity = sample_df.ethnicity.values
    ancestry = sample_df.ancestry.values
    
    ethnicity_names = sample_df.ethnicity.unique()
    eth_le = preprocessing.LabelEncoder().fit(ethnicity_names)
    
    eth_labels = eth_le.classes_
    anc_labels = anc_le.classes_

    plt.figure(figsize=(6,5))
    ax = plt.subplot(111)
    C = confusion_matrix(anc_le.transform(ancestry),
                         eth_le.transform(ethnicity))
    sns.heatmap(C, annot=True, fmt='d', ax=ax)
    
    ax.xaxis.set_ticklabels(eth_labels)
    ax.yaxis.set_ticklabels(anc_labels)
    
    ax.set_xlabel("Self-reported ethnicity")
    ax.set_ylabel("Predicted genetic ancestry")
    ax.collections[0].colorbar.set_label("Num. of individuals")
    plt.title("Ethnicity vs. Ancestry confusion matrix")
    plt.tight_layout()
    plt.savefig(output_filename, dpi=300)
    return

def main(args):
    clf, le, = load_classifier(args.classifier, args.labeler)
    G, sample_df = encode_test_data(args.input, args.aim_dataframe)
    
    pred = le.inverse_transform(clf.predict(G))
    sample_df['ancestry'] = pred
    
    if args.make_confusion_matrix:
        cm_output = args.input.replace(".csv", "") + "_confusion.png"
        test_le = test_confusion_matrix(sample_df, le, cm_output)

    output_df = sample_df[['sample', 'ancestry']]
    output_df.to_csv(args.output, index=False)
    return

if __name__=="__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-c", "--classifier",
        help="Path to trained classifier model",
        default="data/ancestry_svc_clf.pkl")
    
    parser.add_argument("-l", "--labeler",
        help="Path to trained classifier model",
        default="data/ancestry_label_encoder.pkl")
    
    parser.add_argument("-a", "--aim_dataframe",
        help="Path to 1KG AIM dataframe",
        default="data/aim_variants.pkl")
    
    parser.add_argument("-i", "--input",
        help="Path to CSV file containing test individuals")

    parser.add_argument("-o", "--output",
        help="Path to CSV file containing sample names and predicted ancestries",
        default="predicted_ancestries.csv")

    parser.add_argument("-cm", "--make_confusion_matrix",
        help="Use self-reported ethnicity in test samples to make confusion matrix",
        action='store_true')

    main(parser.parse_args())
