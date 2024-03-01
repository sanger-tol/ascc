#!/usr/bin/env python3

"""
Script for kmer count dimensionality reduction, using multiple methods.

Developed by Eerik Aunin (ea10@sanger.ac.uk)

Modified by Damon-Lee Pointon (dp24)
"""

# https://machinelearningmastery.com/dimensionality-reduction-algorithms-with-python/
# https://scikit-learn.org/stable/auto_examples/manifold/plot_compare_methods.html#sphx-glr-auto-examples-manifold-plot-compare-methods-py

import os
matplot_config_path = os.getcwd() + "/matplotconfig/"
fonts_config_path = os.getcwd() + "/fontconfig/"
numba_config_path = os.getcwd() + "/numbaconfig/"

for i in [fonts_config_path, matplot_config_path, numba_config_path]:
    if not os.path.exists(i):
        os.makedirs(i)

os.environ['OSFONTDIR'] = fonts_config_path
os.environ['MPLCONFIGDIR'] = matplot_config_path
os.environ['NUMBA_CACHE_DIR'] = numba_config_path

os.chmod(fonts_config_path, 0o0755)
os.chmod(matplot_config_path, 0o0755)
os.chmod(numba_config_path, 0o0755)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from sklearn import preprocessing
import pandas as pd

import umap
from sklearn import decomposition
from sklearn import manifold
from sklearn import ensemble

import tensorflow as tf
import random

import sys
import argparse
from pathlib import Path
from datetime import datetime
from functools import reduce





# plt.style.use("ggplot")


def reset_random_seeds():
    # https://stackoverflow.com/questions/32419510/how-to-get-reproducible-results-in-keras
    os.environ["PYTHONHASHSEED"] = str(1)
    tf.random.set_seed(1)
    np.random.seed(1)
    random.seed(1)


def print_timestamp(process_name):
    sys.stderr.write("{}: running {}\n".format(datetime.now(), process_name))


def embedding_to_dataframe(embedding, seq_names, plot_title):
    """
    Function for taking a dimensionality reduction embedding and converting it into a pandas dataframe
    """
    plot_title_underscores = plot_title.replace(" ", "_")
    embedding_df = pd.DataFrame(embedding)
    embedding_df.columns = ["embedding_x_" + plot_title_underscores, "embedding_y_" + plot_title_underscores]
    embedding_df["scaff"] = seq_names
    return embedding_df


def load_data(kmer_counts_file):
    """
    Loads kmer counts from an input file, returns the counts as a pandas dataframe and sequence names as a list
    """
    if os.path.isfile(kmer_counts_file) is False:
        sys.stderr.write("kmer counts file ({}) was not found\n".format(kmer_counts_file))
        sys.exit(1)
    if os.stat(kmer_counts_file).st_size == 0:
        sys.stderr.write("kmer counts file ({}) is empty\n".format(kmer_counts_file))
        sys.exit(1)
    data = pd.read_csv(kmer_counts_file, sep=",", index_col=0)
    df_row_count = data.shape[0]
    if df_row_count == 0:
        sys.stderr.write("The kmer counts table ({}) has no rows\n".format(kmer_counts_file))
        sys.exit(1)
    data_to_cluster = data.drop(["seq_len"], axis=1)
    df = data_to_cluster

    seq_names = df.index.values.tolist()
    scalar = preprocessing.MinMaxScaler()
    scalar.fit(data_to_cluster)
    df = scalar.transform(data_to_cluster)
    return df, seq_names


def make_loss_plot(model_history, out_folder, activation_mode):
    """
    Creates the loss vs epoch plot and saves it as a file
    """
    plt.figure(figsize=(12, 6))
    plt.plot(model_history.history["loss"])
    plt.title("Autoencoder {}: Loss vs. Epoch".format(activation_mode))
    plt.ylabel("Loss")
    plt.xlabel("Epoch")
    plt.grid(True)
    plt.savefig(out_folder + "/autoencoder_{}_loss_vs_epoch_plot.png".format(activation_mode))
    plt.close()


def run_autoencoder(df, nr_of_epochs, activation_mode, out_folder):
    """
    Function for running autoencoder for dimensionality reduction of kmer counts
    """
    input_dim = df.shape[1]
    # This is the dimension of the latent space (encoding space)
    latent_dim = 2

    encoder = tf.keras.models.Sequential(
        [
            tf.keras.layers.Dense(128, activation=activation_mode, input_shape=(input_dim,)),
            tf.keras.layers.Dense(64, activation=activation_mode),
            tf.keras.layers.Dense(32, activation=activation_mode),
            tf.keras.layers.Dense(latent_dim, activation=activation_mode),
        ]
    )

    decoder = tf.keras.models.Sequential(
        [
            tf.keras.layers.Dense(64, activation=activation_mode, input_shape=(latent_dim,)),
            tf.keras.layers.Dense(128, activation=activation_mode),
            tf.keras.layers.Dense(256, activation=activation_mode),
            tf.keras.layers.Dense(input_dim, activation=None),
        ]
    )

    autoencoder = tf.keras.models.Model(inputs=encoder.input, outputs=decoder(encoder.output))
    autoencoder.compile(loss="mse", optimizer="adam")
    model_history = autoencoder.fit(df, df, epochs=nr_of_epochs, batch_size=32, verbose=0)
    if out_folder is not None:
        make_loss_plot(model_history, out_folder, activation_mode)
    embedding = encoder(df)
    embedding = embedding.numpy()
    return embedding


def run_tsne(df):
    """
    Function for running t-SNE for dimensionality reduction of kmer counts
    """
    print_timestamp("t-SNE")
    t_sne_perplexity = 30
    df_row_count = df.shape[0]
    if df_row_count < 2:
        sys.stderr.write("Cannot run t-SNE because the kmer counts dataframe has only {} rows\n".format(df_row_count))
    else:
        if t_sne_perplexity >= df_row_count:
            t_sne_perplexity = round(df_row_count / 2)  # t-SNE perplexity must be less than n_samples
            if t_sne_perplexity < 1:
                t_sne_perplexity = 1
            sys.stderr.write("Set t-SNE perplexity value to {}\n".format(t_sne_perplexity))

    embedding = manifold.TSNE(n_components=2).fit_transform(df)
    return embedding

def run_dim_reduction(df, selected_method, n_neighbors_setting=-1, autoencoder_epochs_count=-1, out_folder=None):
    embedding = None
    embedding_title = None
    if selected_method == "pca":
        embedding_title = "PCA"
        print_timestamp(embedding_title)
        pca = decomposition.PCA(n_components=2)
        embedding = pca.fit_transform(df)

    elif selected_method == "umap":
        embedding_title = "UMAP"
        print_timestamp(embedding_title)
        reducer = umap.UMAP(
            random_state=123, n_neighbors=n_neighbors_setting, min_dist=0.1, n_components=2, metric="euclidean"
        )
        embedding = reducer.fit_transform(df)

    elif selected_method == "t-sne":
        embedding_title = "t-SNE"
        embedding = run_tsne(df)

    elif selected_method == "isomap":
        embedding_title = "isomap"
        print_timestamp(embedding_title)
        embedding = manifold.Isomap(n_neighbors=n_neighbors_setting, n_components=2).fit_transform(df)

    elif selected_method == "lle_standard":
        embedding_title = "LLE standard"
        print_timestamp(embedding_title)
        embedding = manifold.LocallyLinearEmbedding(
            n_neighbors=n_neighbors_setting, n_components=2, eigen_solver="dense", method="standard"
        ).fit_transform(df)

    elif selected_method == "lle_hessian":
        embedding_title = "LLE standard"
        print_timestamp(embedding_title)
        embedding = manifold.LocallyLinearEmbedding(
            n_neighbors=n_neighbors_setting, n_components=2, eigen_solver="dense", method="hessian"
        ).fit_transform(df)

    elif selected_method == "lle_modified":
        embedding_title = "LLE modified"
        print_timestamp(embedding_title)
        embedding = manifold.LocallyLinearEmbedding(
            n_neighbors=n_neighbors_setting, n_components=2, eigen_solver="dense", method="modified"
        ).fit_transform(df)

    elif selected_method == "mds":
        embedding_title = "MDS"
        print_timestamp(embedding_title)
        embedding = manifold.MDS(max_iter=100, n_init=1, n_components=2).fit_transform(df)

    elif selected_method == "se":
        embedding_title = "SE"
        print_timestamp(embedding_title)
        embedding = manifold.SpectralEmbedding(n_neighbors=n_neighbors_setting, n_components=2).fit_transform(df)

    elif selected_method == "random_trees":
        # https://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html#sphx-glr-auto-examples-manifold-plot-lle-digits-py
        embedding_title = "random_trees"
        print_timestamp(embedding_title)
        hasher = ensemble.RandomTreesEmbedding(n_estimators=200, random_state=0, max_depth=5)
        x_transformed = hasher.fit_transform(df)
        pca = decomposition.TruncatedSVD(n_components=2)
        embedding = pca.fit_transform(x_transformed)

    elif selected_method == "kernel_pca":
        # https://www.tutorialspoint.com/scikit_learn/scikit_learn_dimensionality_reduction_using_pca.htm
        embedding_title = "KernelPCA"
        print_timestamp(embedding_title)
        kernel_pca = decomposition.KernelPCA(n_components=2, kernel="sigmoid")
        embedding = kernel_pca.fit_transform(df)

    elif selected_method == "pca_svd":
        # https://www.tutorialspoint.com/scikit_learn/scikit_learn_dimensionality_reduction_using_pca.htm
        embedding_title = "PCA with SVD solver"
        print_timestamp(embedding_title)
        pca = decomposition.PCA(n_components=2, svd_solver="randomized")
        embedding = pca.fit_transform(df)

    elif selected_method == "autoencoder_sigmoid":
        print_timestamp("Autoencoder sigmoid")
        # https://ekamperi.github.io/machine%20learning/2021/01/21/encoder-decoder-model.html
        embedding = run_autoencoder(df, autoencoder_epochs_count, "sigmoid", out_folder)

    elif selected_method == "autoencoder_linear":
        # https://ekamperi.github.io/machine%20learning/2021/01/21/encoder-decoder-model.html
        embedding_title = "Autoencoder linear"
        print_timestamp(embedding_title)
        embedding = run_autoencoder(df, autoencoder_epochs_count, "linear", out_folder)

    elif selected_method == "autoencoder_tanh":
        # https://ekamperi.github.io/machine%20learning/2021/01/21/encoder-decoder-model.html
        embedding_title = "Autoencoder tanh"
        print_timestamp(embedding_title)
        embedding = run_autoencoder(df, autoencoder_epochs_count, "tanh", out_folder)

    elif selected_method == "autoencoder_selu":
        # https://ekamperi.github.io/machine%20learning/2021/01/21/encoder-decoder-model.html
        embedding_title = "Autoencoder selu"
        print_timestamp(embedding_title)
        embedding = run_autoencoder(df, autoencoder_epochs_count, "selu", out_folder)

    elif selected_method == "autoencoder_relu":
        # https://ekamperi.github.io/machine%20learning/2021/01/21/encoder-decoder-model.html
        embedding_title = "Autoencoder relu"
        print_timestamp(embedding_title)
        embedding = run_autoencoder(df, autoencoder_epochs_count, "relu", out_folder)

    elif selected_method == "nmf":
        # https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html
        embedding_title = "Non-Negative Matrix Factorization"
        print_timestamp(embedding_title)
        nmf_model = decomposition.NMF(n_components=2, init="random", random_state=0)
        embedding = nmf_model.fit_transform(df)
    return embedding, embedding_title


def main(kmer_counts_file, out_folder, selected_methods, n_neighbors_setting, autoencoder_epochs_count):
    reset_random_seeds()

    df, seq_names = load_data(kmer_counts_file)
    df_row_count = df.shape[0]
    if df_row_count == 1:
        sys.stderr.write(
            "Skipping the dimensionality reduction of kmer counts, as the kmer counts table has only one row"
        )
        # Generate an empty file to satisfy nextflow expecting a file from script finishing with no file with small output
        with open("EMPTY_kmers_dim_reduction_embeddings.csv") as empty_file:
            empty_file.write("FILE TO SMALL FOR ANALYSIS")
        sys.exit(0)

    Path(out_folder).mkdir(parents=True, exist_ok=True)
    selected_methods = selected_methods.split(",")

    if autoencoder_epochs_count == -1:
        autoencoder_epochs_count = df_row_count

    embeddings_list = list()
    for selected_method in selected_methods:
        embedding, embedding_title = run_dim_reduction(
            df,
            selected_method,
            n_neighbors_setting=n_neighbors_setting,
            autoencoder_epochs_count=autoencoder_epochs_count,
            out_folder=out_folder,
        )
        embedding_df = embedding_to_dataframe(embedding, seq_names, embedding_title)
        embeddings_list.append(embedding_df)

    out_df = reduce(
        lambda left, right: pd.merge(left, right, on=["scaff"], how="outer"),  # Merge DataFrames in list
        embeddings_list,
    )
    out_path = out_folder + "/kmers_dim_reduction_embeddings.csv"
    out_df.to_csv(out_path, index=False)

    os._exit(0)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--version", action="version", version="1.0")
    parser.add_argument("kmer_counts_file", type=str, help="Path to input CSV file with kmer counts")
    parser.add_argument("out_folder", type=str, help="Path to folder where output files will be written")
    parser.add_argument(
        "--selected_methods",
        type=str,
        help="Comma separated string with the selected dimensionality reduction methods",
        default="pca",
    )
    parser.add_argument(
        "--n_neighbors_setting",
        type=int,
        help="n_neighbors parameter value for the methods that have this parameter (default: 13)",
        default=13,
    )
    parser.add_argument(
        "--autoencoder_epochs_count",
        type=int,
        help="Autoencoder epochs count (default: assign automatically as 3x number of sequences in input)",
        default=-1,
    )
    args = parser.parse_args()
    main(
        args.kmer_counts_file,
        args.out_folder,
        args.selected_methods,
        args.n_neighbors_setting,
        args.autoencoder_epochs_count,
    )
