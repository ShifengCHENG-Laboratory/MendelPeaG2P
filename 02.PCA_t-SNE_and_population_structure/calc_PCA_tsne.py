import pandas as pd
import numpy as np
from sklearn.preprocessing import OneHotEncoder
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import sys
import os

def calc_pca(df, output_prefix):
    pca = PCA(n_components=3)
    data_pca = pca.fit_transform(df)

    # create a dataframe for PCA results
    df_pca = pd.DataFrame(data=data_pca, columns=[f"PC{i+1}" for i in range(3)], index=df.index)
    df_pca.reset_index(inplace=True)
    df_pca.rename(columns={'index': 'Sample'}, inplace=True)

    # save PCA results to a csv file
    df_pca.to_csv(f"{output_prefix}_PCA_results.csv", index=False)

    # create a dataframe for the eigenvalues (explained variance) and explained variance ratio (percentage of variance explained by each of the selected components)
    df_eigenvalues_and_ratio = pd.DataFrame({
        'PCs': [f"PC{i+1}" for i in range(3)],
        'Eigenvalue': pca.explained_variance_,
        'Variance Percent (%)': pca.explained_variance_ratio_ * 100
    })

    # output to a file
    df_eigenvalues_and_ratio.to_csv(f"{output_prefix}_PCA_eigenvalues.tsv", sep="\t", index=False)

def calc_tsne(df, output_prefix):
    # perform T-SNE
    tsne = TSNE(n_components=2, init='pca', random_state=0)
    tsne_results = tsne.fit_transform(df)

    # create a dataframe for T-SNE results
    df_tsne = pd.DataFrame(data=tsne_results, columns=["Dim1", "Dim2"], index=df.index)
    df_tsne.reset_index(inplace=True)
    df_tsne.rename(columns={'index': 'Sample'}, inplace=True)

    # save T-SNE results to a csv file
    df_tsne.to_csv(f"{output_prefix}_tSNE_results.csv", index=False)

def main(input_matrix, sample_list_file):
    # load your matrix
    df = pd.read_csv(input_matrix, sep="\t", index_col=0)

    # load sample list
    with open(sample_list_file, 'r') as f:
        sample_list = f.read().splitlines()

    # filter df based on sample list
    df = df.loc[:, df.columns.intersection(sample_list)]

    # transpose the dataframe
    df = df.transpose()

    output_prefix = os.path.splitext(input_matrix)[0]

    # without one-hot encoding
    calc_pca(df, output_prefix + "_no_one_hot")
    calc_tsne(df, output_prefix + "_no_one_hot")

    # with one-hot encoding
    enc = OneHotEncoder(sparse=False)
    df_onehot = pd.DataFrame(enc.fit_transform(df), columns=enc.get_feature_names(df.columns), index=df.index)
    
    calc_pca(df_onehot, output_prefix + "_one_hot")
    calc_tsne(df_onehot, output_prefix + "_one_hot")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python3 calc_pca_and_tsne.py in.matrix.txt sample_list.txt")
        sys.exit(1)

    main(sys.argv[1], sys.argv[2])

