#!/usr/bin/env python3

"""
This script performs dimensionality reduction on k-mer count tables using various methods.

Developed by Eerik Aunin (eeaunin@gmail.com)

Modified by Damon-Lee Pointon (dp24)

The autoencoder part is partially based on:
    "Disentangling cobionts and contamination in long-read genomic data using sequence composition",
        Claudia C Weber,
        G3 Genes|Genomes|Genetics,
        2024.

The preprocessing pipeline includes:
1. Normalisation by sequence length (using the 'seq_len' column if present).
2. Scaling based on the selected method:
    - PCA, LLE, Kernel PCA: Standard scaling (mean=0, std=1).
    - UMAP, t-SNE, Isomap: Log scaling for skewed distributions.
    - Random Trees, NMF: No scaling (raw data is used).
    - Autoencoder: tf-idf scaling

"""

# https://machinelearningmastery.com/dimensionality-reduction-algorithms-with-python/
# https://scikit-learn.org/stable/auto_examples/manifold/plot_compare_methods.html#sphx-glr-auto-examples-manifold-plot-compare-methods-py

import os

matplot_config_path = os.path.join(os.getcwd(), "matplotconfig")
fonts_config_path = os.path.join(os.getcwd(), "fontconfig")
numba_config_path = os.path.join(os.getcwd(), "numbaconfig")

for i in [fonts_config_path, matplot_config_path, numba_config_path]:
    if not os.path.exists(i):
        os.makedirs(i, exist_ok=True)

os.environ["OSFONTDIR"] = fonts_config_path
os.environ["MPLCONFIGDIR"] = matplot_config_path
os.environ["NUMBA_CACHE_DIR"] = numba_config_path

os.chmod(fonts_config_path, 0o0775)
os.chmod(matplot_config_path, 0o0775)
os.chmod(numba_config_path, 0o0775)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

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
import gc  # for garbage collection
from sklearn.preprocessing import MinMaxScaler, StandardScaler
from sklearn.neighbors import NearestNeighbors
from tensorflow.keras import backend as K
from sklearn.metrics import (
    silhouette_score,
    calinski_harabasz_score,
    davies_bouldin_score,
)
from sklearn.cluster import HDBSCAN
from sklearn.metrics import adjusted_rand_score


def reset_random_seeds():
    # https://stackoverflow.com/questions/32419510/how-to-get-reproducible-results-in-keras
    os.environ["PYTHONHASHSEED"] = str(1)
    tf.random.set_seed(1)
    np.random.seed(1)
    random.seed(1)


def print_timestamp(process_name):
    sys.stderr.write("{}: running {}\n".format(datetime.now(), process_name))


def log_message(message: str, level: str = "warning", timestamp: bool = True) -> None:
    """
    Enhanced logging function with timestamps and consistent formatting.
    """
    levels = {
        "info": "[INFO]",
        "warning": "[WARNING]",
        "error": "[ERROR]",
        "debug": "[DEBUG]",
    }
    prefix = levels.get(level.lower(), "[LOG]")
    timestamp_str = (
        f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] " if timestamp else ""
    )
    sys.stderr.write(f"{timestamp_str}{prefix} {message}\n")


def embedding_to_dataframe(embedding, seq_names, plot_title):
    """
    Function for taking a dimensionality reduction embedding and converting it into a pandas dataframe
    """
    plot_title_underscores = plot_title.replace(" ", "_")
    embedding_df = pd.DataFrame(embedding)

    # Create column names for all dimensions
    dim_names = []
    for i in range(embedding.shape[1]):
        dim_names.append(f"embedding_dim_{i+1}_{plot_title_underscores}")

    embedding_df.columns = dim_names
    embedding_df["scaff"] = seq_names
    return embedding_df


def load_data(kmer_counts_file):
    """
    Loads k-mer counts from an input file, returns the counts as a pandas DataFrame,
    sequence names as a list, and sequence lengths as a Series.
    """

    if os.path.isfile(kmer_counts_file) is False:
        sys.stderr.write(
            "kmer counts file ({}) was not found\n".format(kmer_counts_file)
        )
        sys.exit(1)
    if os.stat(kmer_counts_file).st_size == 0:
        sys.stderr.write("kmer counts file ({}) is empty\n".format(kmer_counts_file))
        sys.exit(1)

    try:
        data = pd.read_csv(kmer_counts_file, sep=",", index_col=0)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        sys.stderr.write("The input file is empty or invalid\n")
        sys.exit(1)

    if data.isnull().values.any():
        sys.stderr.write(
            "Warning: Input data contains NaN values. These will be handled by the preprocessing step\n"
        )

    try:
        seq_lengths = data["seq_len"]
        data_to_cluster = data.drop(["seq_len"], axis=1)
    except KeyError:
        sys.stderr.write("Warning: 'seq_len' column not found in input data\n")
        seq_lengths = pd.Series(1, index=data.index)  # Default to 1 if missing
        data_to_cluster = data

    if seq_lengths.eq(1).sum() > 0:
        log_message(
            f"Warning: {seq_lengths.eq(1).sum()} sequences are missing 'seq_len' values. Defaulting to 1.",
            level="warning",
        )

    seq_names = data_to_cluster.index.values.tolist()

    # Dimension check
    if data_to_cluster.shape[0] != len(seq_names):
        log_message(
            f"Mismatch between number of sequences ({len(seq_names)}) and data rows ({data_to_cluster.shape[0]})",
            level="error",
        )
        sys.exit(1)

    return data_to_cluster, seq_names, seq_lengths


def remove_constant_features(df, variance_threshold=1e-5):
    """
    Remove features that are constant or near-constant
    """
    # Calculate variances
    variances = df.var()
    near_constant = variances <= variance_threshold
    n_removed = sum(near_constant)

    if n_removed > 0:
        percent_removed = (n_removed / df.shape[1]) * 100
        log_message(
            f"Removing {n_removed} near-constant k-mers "
            f"({percent_removed:.2f}% of features)",
            level="info",
        )
        df = df.loc[:, ~near_constant]

    return df, n_removed


def calculate_tf_idf(df):
    """
    Calculate TF-IDF for k-mer counts with numerical stability safeguards.

    Parameters:
    df: DataFrame with k-mer counts (rows are sequences, columns are k-mers)

    Returns:
    DataFrame with TF-IDF values
    """
    # Calculate TF (term frequency)
    # Add small epsilon to prevent division by zero
    row_sums = df.sum(axis=1)
    epsilon = 1e-10
    row_sums = row_sums.clip(lower=epsilon)  # Ensure no zeros in denominators
    tf = df.div(row_sums, axis=0)

    # Calculate IDF (inverse document frequency)
    n_sequences = len(df)
    # Add epsilon to both numerator and denominator for numerical stability
    document_frequency = (df > 0).sum(axis=0) + epsilon
    idf = np.log((n_sequences + epsilon) / document_frequency)

    # Calculate TF-IDF
    tf_idf = tf.multiply(idf)

    # Final cleanup of any remaining numerical instabilities
    tf_idf = tf_idf.clip(lower=0.0)  # Ensure no negative values

    # Replace any remaining NaN or inf values with 0
    tf_idf = tf_idf.replace([np.inf, -np.inf], 0.0)
    tf_idf = tf_idf.fillna(0.0)
    return tf_idf


def preprocess_data(df, seq_lengths, scaling_method="log", use_tf_idf=False):
    data_array = None
    if use_tf_idf:
        # TF-IDF inherently handles length normalisation
        data_array = calculate_tf_idf(df).values

        # Additional numerical stability
        data_array = np.nan_to_num(data_array, nan=0.0, posinf=0.0, neginf=0.0)

        # Scale to avoid extreme values
        if np.any(data_array):  # Only if we have non-zero values
            scaler = MinMaxScaler(feature_range=(1e-6, 1.0 - 1e-6))
            data_array = scaler.fit_transform(data_array)

        return data_array
    else:
        # For non-TF-IDF approaches
        df_normalised = df.div(
            seq_lengths.clip(lower=1e-10), axis=0
        )  # Prevent division by zero
        data_array = df_normalised.values

        # Apply scaling
        if scaling_method == "log":
            return np.log1p(data_array)
        elif scaling_method == "minmax":
            scaler = MinMaxScaler(feature_range=(1e-6, 1.0 - 1e-6))
            return scaler.fit_transform(data_array)
        elif scaling_method == "standard":
            scaler = StandardScaler()
            return scaler.fit_transform(data_array)
        else:
            # No further scaling
            return data_array


def transform_kmer_counts(df):
    """
    Transform k-mer counts to handle varying magnitudes via log(1 + x).
    """
    # Add small constant before log to handle zeros
    return np.log1p(df)


def validate_n_neighbors(n_neighbors, n_samples, method_name):
    """
    Validates and adjusts n_neighbors parameter for methods that use it.
    Returns adjusted n_neighbors value.
    """
    if n_neighbors >= n_samples:
        adjusted = n_samples - 1
        log_message(
            f"Reducing n_neighbors from {n_neighbors} to {adjusted} for {method_name} "
            f"due to dataset size ({n_samples} samples)",
            level="warning",
        )
        return adjusted
    return n_neighbors


def evaluate_clustering(embedding, labels):
    """
    Evaluate clustering quality using multiple metrics
    """
    metrics = {}

    n_samples = len(labels)
    n_unique_labels = len(set(labels))

    # Check if we have enough samples and labels for meaningful metrics
    if n_samples < 4 or n_unique_labels < 2 or n_unique_labels >= n_samples:
        log_message(
            f"Dataset too small for clustering metrics (samples: {n_samples}, "
            f"unique labels: {n_unique_labels}). Skipping metrics calculation.",
            level="warning",
        )
        return {
            "silhouette": float("nan"),
            "calinski_harabasz": float("nan"),
            "davies_bouldin": float("nan"),
            "adjusted_rand": float("nan"),
        }

    # Silhouette score (-1 to 1, higher is better)
    try:
        metrics["silhouette"] = silhouette_score(embedding, labels)
    except Exception as e:
        log_message(f"Failed to calculate silhouette score: {str(e)}", level="warning")
        metrics["silhouette"] = float("nan")

    # Calinski-Harabasz Index (higher is better)
    try:
        metrics["calinski_harabasz"] = calinski_harabasz_score(embedding, labels)
    except Exception as e:
        log_message(
            f"Failed to calculate Calinski-Harabasz score: {str(e)}", level="warning"
        )
        metrics["calinski_harabasz"] = float("nan")

    # Davies-Bouldin Index (lower is better)
    try:
        metrics["davies_bouldin"] = davies_bouldin_score(embedding, labels)
    except Exception as e:
        log_message(
            f"Failed to calculate Davies-Bouldin score: {str(e)}", level="warning"
        )
        metrics["davies_bouldin"] = float("nan")

    # HDBSCAN clustering comparison
    try:
        clusterer = HDBSCAN(min_cluster_size=min(5, n_samples - 1))
        cluster_labels = clusterer.fit_predict(embedding)
        metrics["adjusted_rand"] = adjusted_rand_score(labels, cluster_labels)
    except Exception as e:
        log_message(
            f"Failed to calculate adjusted rand index: {str(e)}", level="warning"
        )
        metrics["adjusted_rand"] = float("nan")

    return metrics


def find_best_n_neighbors_setting(
    data,
    method_name,
    n_neighbors_range=[5, 10, 15, 20, 30],
    min_cluster_sizes=[5, 10, 15],
    max_optimisation_samples=1000,
):
    """
    Optimise n_neighbors parameter for different dimensionality reduction methods.

    Parameters:
        data: array-like, input data
        method_name: str, one of ['umap', 'isomap', 'lle_standard', 'lle_hessian', 'lle_modified', 'se']
        n_neighbors_range: list of int, n_neighbors values to try
        min_cluster_sizes: list of int, HDBSCAN min_cluster_size values to try
        max_optimisation_samples: int, maximum number of samples to use for optimisation

    Returns:
        tuple: (best_n_neighbors, optimisation_results_dict)
    """

    # Add validation
    if len(data) < min(n_neighbors_range):
        log_message(
            f"Dataset too small for n_neighbors optimisation. Using minimum value.",
            level="warning",
        )
        return min(n_neighbors_range), {}

    # Sample data if needed
    if len(data) > max_optimisation_samples:
        indices = np.random.choice(len(data), max_optimisation_samples, replace=False)
        optimisation_data = data[indices]
    else:
        optimisation_data = data

    best_score = -1
    best_n_neighbors = None
    results = {}

    for n_neighbors in n_neighbors_range:
        # Apply dimensionality reduction based on method
        try:
            if method_name == "umap":
                embedding = umap.UMAP(n_neighbors=n_neighbors).fit_transform(
                    optimisation_data
                )
            elif method_name == "isomap":
                embedding = manifold.Isomap(n_neighbors=n_neighbors).fit_transform(
                    optimisation_data
                )
            elif method_name == "lle_standard":
                embedding = manifold.LocallyLinearEmbedding(
                    n_neighbors=n_neighbors, method="standard"
                ).fit_transform(optimisation_data)
            elif method_name == "lle_hessian":
                embedding = manifold.LocallyLinearEmbedding(
                    n_neighbors=n_neighbors, method="hessian"
                ).fit_transform(optimisation_data)
            elif method_name == "lle_modified":
                embedding = manifold.LocallyLinearEmbedding(
                    n_neighbors=n_neighbors, method="modified"
                ).fit_transform(optimisation_data)
            elif method_name == "se":
                embedding = manifold.SpectralEmbedding(
                    n_neighbors=n_neighbors
                ).fit_transform(optimisation_data)
            else:
                raise ValueError(f"Unknown method: {method_name}")

            # Try different HDBSCAN parameters
            method_scores = []
            for min_cluster_size in min_cluster_sizes:
                clusterer = HDBSCAN(min_cluster_size=min_cluster_size)
                cluster_labels = clusterer.fit_predict(embedding)

                if len(np.unique(cluster_labels)) > 1:
                    score = silhouette_score(embedding, cluster_labels)
                    method_scores.append(score)

            # Average score across different cluster sizes
            if method_scores:
                avg_score = np.mean(method_scores)
                results[n_neighbors] = avg_score

                if avg_score > best_score:
                    best_score = avg_score
                    best_n_neighbors = n_neighbors

        except Exception as e:
            log_message(
                f"Error optimising {method_name} with n_neighbors={n_neighbors}: {str(e)}",
                level="warning",
            )
            continue

    if best_n_neighbors is None:
        # Provide default values based on method
        n_components = min(data.shape[1], 5)  # assuming default n_components=5

        if method_name.startswith("lle_"):
            lle_method = method_name.split("_")[
                1
            ]  # extract 'standard', 'hessian', or 'modified'
            min_neighbors = get_min_neighbors_for_lle(lle_method, n_components)
            best_n_neighbors = max(
                min_neighbors, 15
            )  # Use 15 or minimum required, whichever is larger
            log_message(
                f"n_neighbors ptimisation failed for {method_name}. Using default n_neighbors={best_n_neighbors} "
                f"(minimum required for this method: {min_neighbors})",
                level="warning",
            )
        else:
            best_n_neighbors = 15  # Default fallback for other methods
            log_message(
                f"n_neighbors optimisation failed for {method_name}. Using default n_neighbors=15",
                level="warning",
            )
    return best_n_neighbors, results


def check_memory_requirements(df, selected_methods):
    """
    Estimate memory requirements for different methods
    """
    n_samples, n_features = df.shape

    # Different methods have different memory scaling characteristics
    memory_warnings = {
        # Manifold learning methods with quadratic memory scaling
        "t-sne": (
            10000,
            "quadratic memory scaling with the number of samples. Try reducing the number of samples or using a lower perplexity.",
        ),
        "mds": (
            5000,
            "quadratic memory scaling with the number of samples. Consider reducing dataset size.",
        ),
        "isomap": (
            10000,
            "quadratic memory scaling due to geodesic distance calculations. Consider reducing n_neighbors or dataset size.",
        ),
        "lle_standard": (
            10000,
            "quadratic memory scaling with the number of samples. Consider reducing n_neighbors.",
        ),
        "lle_hessian": (
            5000,
            "high memory usage due to Hessian calculations. Consider reducing n_neighbors.",
        ),
        "lle_modified": (
            10000,
            "quadratic memory scaling with the number of samples. Consider reducing n_neighbors.",
        ),
        "se": (
            10000,
            "quadratic memory scaling due to eigendecomposition of the similarity matrix.",
        ),
        # Linear methods with better memory scaling
        "pca": (
            50000,
            "linear memory scaling but may require significant memory for very high-dimensional data.",
        ),
        "kernel_pca": (
            20000,
            "memory usage depends on kernel type, may be significant for large datasets.",
        ),
        "pca_svd": (
            50000,
            "randomized solver helps with memory usage but still requires significant memory for large datasets.",
        ),
        "nmf": (
            30000,
            "iterative nature requires storing multiple copies of the data matrix.",
        ),
        # Tree-based method
        "random_trees": (
            100000,
            "memory scales with n_estimators and max_depth parameters.",
        ),
        # UMAP is generally memory-efficient
        "umap": (
            100000,
            "relatively memory-efficient but may require significant memory for very large datasets.",
        ),
        # Autoencoder variants
        "autoencoder_relu": (
            20000,
            "significant GPU/CPU memory required for training. Consider reducing batch size.",
        ),
        "autoencoder_tanh": (
            20000,
            "significant GPU/CPU memory required for training. Consider reducing batch size.",
        ),
        "autoencoder_sigmoid": (
            20000,
            "significant GPU/CPU memory required for training. Consider reducing batch size.",
        ),
        "autoencoder_linear": (
            20000,
            "significant GPU/CPU memory required for training. Consider reducing batch size.",
        ),
        "autoencoder_selu": (
            20000,
            "significant GPU/CPU memory required for training. Consider reducing batch size.",
        ),
    }

    # Additional warning for high-dimensional data
    if n_features > 10000:
        log_message(
            f"Warning: Input data has {n_features} features. "
            "This might require significant memory for all methods.",
            level="warning",
        )

    for method in selected_methods:
        if method in memory_warnings:
            sample_threshold, warning_msg = memory_warnings[method]
            if n_samples > sample_threshold:
                log_message(
                    f"Warning: {method} {warning_msg} "
                    f"Current dataset has {n_samples} samples.",
                    level="warning",
                )


class DataGenerator(tf.keras.utils.Sequence):
    """
    Memory-efficient data generator with dynamic batch sizing. Used for autoencoder
    """

    def __init__(self, data, batch_size=32, shuffle=True, max_memory_gb=4):
        # Validation of input parameters
        if batch_size <= 0:
            raise ValueError("batch_size must be positive")
        if max_memory_gb <= 0:
            raise ValueError("max_memory_gb must be positive")
        if not isinstance(data, np.ndarray):
            raise ValueError("data must be a numpy array")

        self.data = data
        self.shuffle = shuffle
        self.indices = np.arange(len(data))

        # Calculate memory-safe batch size
        data_size_bytes = data.nbytes / len(data)  # bytes per sample
        max_memory_bytes = max_memory_gb * 1024**3
        memory_safe_batch = int(
            max_memory_bytes / (data_size_bytes * 3)
        )  # Factor of 3 for safety
        self.batch_size = min(batch_size, memory_safe_batch)

        if self.batch_size != batch_size:
            log_message(
                f"Adjusted batch size from {batch_size} to {self.batch_size} "
                "to prevent memory issues",
                level="info",
            )

        self.on_epoch_end()

    def __getitem__(self, idx):
        batch_indices = self.indices[
            idx * self.batch_size : (idx + 1) * self.batch_size
        ]
        batch_data = self.data[batch_indices]
        # For VAE, input = target
        return batch_data, batch_data

    def __len__(self):
        return int(np.ceil(len(self.data) / self.batch_size))

    def on_epoch_end(self):
        if self.shuffle:
            np.random.shuffle(self.indices)


def run_autoencoder(
    df,
    seq_lengths,
    nr_of_epochs=100,
    activation_mode="relu",
    out_folder=None,
    n_components=10,
    batch_size=None,
    learning_rate=0.001,
    scaling_method="log",
    return_umap=False,
    seq_names=None,
    optimise_n_neighbors=False,
):
    """
    Function for running a variational autoencoder for dimensionality reduction of k-mer counts.
    Variational Autoencoder (VAE) implementation for k-mer count dimensionality reduction.

    This implementation includes:
    - Proper KL divergence warm-up (β-VAE style)
    - Batch normalisation
    - Non-negative output enforcement (via ReLU on decoder)
    - Early stopping and learning rate scheduling based on validation metrics
    - Optional log transform of input k-mer counts
    """

    class VAE(tf.keras.Model):
        def __init__(
            self, input_dim, n_components, activation=activation_mode, **kwargs
        ):
            super().__init__(**kwargs)

            # Initialise metrics
            self.total_loss_tracker = tf.keras.metrics.Mean(name="loss")
            self.reconstruction_loss_tracker = tf.keras.metrics.Mean(
                name="reconstruction_loss"
            )
            self.kl_loss_tracker = tf.keras.metrics.Mean(name="kl_loss")
            self.val_total_loss_tracker = tf.keras.metrics.Mean(name="val_loss")
            self.val_reconstruction_loss_tracker = tf.keras.metrics.Mean(
                name="val_reconstruction_loss"
            )
            self.val_kl_loss_tracker = tf.keras.metrics.Mean(name="val_kl_loss")

            # Dynamic intermediate dimensions
            intermediate_dim = max(32, input_dim // 8)

            # Encoder
            encoder_inputs = tf.keras.layers.Input(shape=(input_dim,))

            # Initial normalisation
            x = tf.keras.layers.BatchNormalization()(encoder_inputs)

            # First dense block with residual connection
            x1 = tf.keras.layers.Dense(intermediate_dim, activation=activation)(x)
            x = tf.keras.layers.Dense(intermediate_dim)(x)
            x = tf.keras.layers.Add()([x1, x])
            x = tf.keras.layers.Activation(activation)(x)
            x = tf.keras.layers.BatchNormalization()(x)
            x = tf.keras.layers.Dropout(0.1)(x)

            # Second dense block
            x2 = tf.keras.layers.Dense(intermediate_dim // 2, activation=activation)(x)
            x = tf.keras.layers.Dense(intermediate_dim // 2)(x)
            x = tf.keras.layers.Add()([x2, x])
            x = tf.keras.layers.Activation(activation)(x)
            x = tf.keras.layers.BatchNormalization()(x)

            # Final encoder layers
            x = tf.keras.layers.Dense(
                intermediate_dim // 4,
                activation=activation,
                activity_regularizer=tf.keras.regularizers.l1(1e-6),
            )(x)

            # Latent space
            z_mean = tf.keras.layers.Dense(n_components, name="z_mean")(x)
            z_log_var = tf.keras.layers.Dense(n_components, name="z_log_var")(x)

            # Decoder
            decoder_inputs = tf.keras.layers.Input(shape=(n_components,))

            # Initial decoder processing
            x = tf.keras.layers.BatchNormalization()(decoder_inputs)

            # Mirror encoder architecture in reverse
            x = tf.keras.layers.Dense(intermediate_dim // 4, activation=activation)(x)

            # First dense block
            x1 = tf.keras.layers.Dense(intermediate_dim // 2, activation=activation)(x)
            x = tf.keras.layers.Dense(intermediate_dim // 2)(x)
            x = tf.keras.layers.Add()([x1, x])
            x = tf.keras.layers.Activation(activation)(x)
            x = tf.keras.layers.BatchNormalization()(x)
            x = tf.keras.layers.Dropout(0.1)(x)

            # Second dense block
            x2 = tf.keras.layers.Dense(intermediate_dim, activation=activation)(x)
            x = tf.keras.layers.Dense(intermediate_dim)(x)
            x = tf.keras.layers.Add()([x2, x])
            x = tf.keras.layers.Activation(activation)(x)
            x = tf.keras.layers.BatchNormalization()(x)

            # Final decoder output
            decoder_outputs = tf.keras.layers.Dense(input_dim, activation="sigmoid")(x)

            # Create encoder and decoder models
            self.encoder = tf.keras.Model(
                encoder_inputs, [z_mean, z_log_var], name="encoder"
            )
            self.decoder = tf.keras.Model(
                decoder_inputs, decoder_outputs, name="decoder"
            )

            # Initialise beta (KL weight) with a smaller value for better training stability
            self.beta = tf.Variable(0.001, trainable=False, dtype=tf.float32)

        @property
        def metrics(self):
            return [
                self.total_loss_tracker,
                self.reconstruction_loss_tracker,
                self.kl_loss_tracker,
                self.val_total_loss_tracker,
                self.val_reconstruction_loss_tracker,
                self.val_kl_loss_tracker,
            ]

        def compute_loss(self, x, training=False):
            z_mean, z_log_var = self.encoder(x)

            batch = tf.shape(z_mean)[0]
            dim = tf.shape(z_mean)[1]
            epsilon = tf.random.normal(shape=(batch, dim))
            z = z_mean + tf.exp(0.5 * z_log_var) * epsilon

            reconstruction = self.decoder(z)
            reconstruction_loss = tf.reduce_mean(tf.square(x - reconstruction))
            reconstruction_loss *= tf.cast(
                tf.shape(x)[1], tf.float32
            )  # Scale by input dimension

            kl_loss = -0.5 * tf.reduce_mean(
                1 + z_log_var - tf.square(z_mean) - tf.exp(z_log_var)
            )

            total_loss = reconstruction_loss + self.beta * kl_loss

            return total_loss, reconstruction_loss, kl_loss

        def train_step(self, data):
            if isinstance(data, tuple):
                x = data[0]
            else:
                x = data

            # Reset metrics
            self.total_loss_tracker.reset_states()
            self.reconstruction_loss_tracker.reset_states()
            self.kl_loss_tracker.reset_states()

            with tf.GradientTape() as tape:
                total_loss, reconstruction_loss, kl_loss = self.compute_loss(
                    x, training=True
                )

            # Compute and clip gradients
            trainable_vars = self.trainable_variables
            gradients = tape.gradient(total_loss, trainable_vars)
            gradients = [
                tf.clip_by_norm(g, 1.0) if g is not None else g for g in gradients
            ]

            # Update weights
            self.optimizer.apply_gradients(zip(gradients, trainable_vars))

            # Update metrics
            self.total_loss_tracker.update_state(total_loss)
            self.reconstruction_loss_tracker.update_state(reconstruction_loss)
            self.kl_loss_tracker.update_state(kl_loss)

            return {
                "loss": self.total_loss_tracker.result(),
                "reconstruction_loss": self.reconstruction_loss_tracker.result(),
                "kl_loss": self.kl_loss_tracker.result(),
            }

        def test_step(self, data):
            if isinstance(data, tuple):
                x = data[0]
            else:
                x = data

            # Reset validation metrics
            self.val_total_loss_tracker.reset_states()
            self.val_reconstruction_loss_tracker.reset_states()
            self.val_kl_loss_tracker.reset_states()

            total_loss, reconstruction_loss, kl_loss = self.compute_loss(
                x, training=False
            )

            # Update validation metrics
            self.val_total_loss_tracker.update_state(total_loss)
            self.val_reconstruction_loss_tracker.update_state(reconstruction_loss)
            self.val_kl_loss_tracker.update_state(kl_loss)

            return {
                "loss": self.val_total_loss_tracker.result(),
                "reconstruction_loss": self.val_reconstruction_loss_tracker.result(),
                "kl_loss": self.val_kl_loss_tracker.result(),
            }

    def create_validation_split(data_array, validation_split=0.2, min_val_samples=50):
        """
        Create a more robust validation split that ensures:
        1. Minimum number of validation samples
        2. Stratification by sequence length (if possible)
        3. Random but reproducible splits
        """
        dataset_size = data_array.shape[0]

        # Adjust validation_split if dataset is small
        if dataset_size < 250:
            min_val_samples = min(min_val_samples, dataset_size // 5)
            validation_split = max(min_val_samples / dataset_size, 0.1)
            log_message(
                f"Small dataset detected ({dataset_size} samples). "
                f"Adjusted validation split to {validation_split:.2f}",
                level="info",
            )

        # Calculate split sizes
        val_size = int(dataset_size * validation_split)
        if val_size < min_val_samples:
            val_size = min(min_val_samples, dataset_size // 2)
            log_message(
                f"Adjusted validation size to {val_size} samples to ensure reliable validation",
                level="info",
            )

        # Create indices for split
        indices = np.random.permutation(dataset_size)
        val_idx = indices[:val_size]
        train_idx = indices[val_size:]

        return train_idx, val_idx

    def determine_batch_size(data_array, min_batch=16, max_batch=256):
        """
        Determine optimal batch size based on dataset characteristics.
        """
        n_samples = len(data_array)

        # Base batch size on dataset size
        if n_samples < 100:
            batch_size = min_batch
        elif n_samples < 1000:
            batch_size = min(64, n_samples // 10)
        elif n_samples < 10000:
            batch_size = min(128, n_samples // 20)
        else:
            batch_size = min(max_batch, n_samples // 50)

        # Adjust for memory constraints
        data_size_gb = data_array.nbytes / (1024**3)
        if data_size_gb > 1:  # If dataset is larger than 1GB
            batch_size = min(batch_size, int(128 / data_size_gb))

        return max(min_batch, batch_size)

    def check_model_architecture(input_dim, n_components):
        """
        Verify model architecture
        """
        if input_dim <= 0 or n_components <= 0:
            raise ValueError("Input and latent dimensions must be positive")
        bottleneck_ratio = input_dim / n_components
        if bottleneck_ratio < 2:
            log_message(
                f"Warning: Bottleneck ratio ({bottleneck_ratio:.1f}) is less than 2. "
                "Consider increasing compression ratio.",
                level="warning",
            )
        if n_components < 2:
            log_message(
                "Warning: Latent dimension less than 2 may be too restrictive.",
                level="warning",
            )

    # -------------------------
    # Before starting: GPU availability check
    # -------------------------
    try:
        # Check if GPU is available
        physical_devices = tf.config.list_physical_devices("GPU")
        if len(physical_devices) == 0:
            log_message("No GPU found. Running on CPU may be slower.", level="warning")
    except:
        log_message("Could not check GPU availability", level="warning")

    # -------------------------
    # 1. Scaling and setting up a generator
    # -------------------------

    data_array = preprocess_data(
        df, seq_lengths, scaling_method=scaling_method, use_tf_idf=True
    )
    if np.any(np.isnan(data_array)) or np.any(np.isinf(data_array)):
        raise ValueError("Input data contains NaN or infinite values.")

    # Determine whether to use generator based on dataset size
    # You can adjust this threshold based on your system's memory

    # Determine optimal batch size if not provided
    if batch_size is None:
        batch_size = determine_batch_size(data_array)
        log_message(f"Using automatically determined batch size: {batch_size}")

    # Use the improved validation split
    train_idx, val_idx = create_validation_split(data_array)

    # Memory-efficient data handling
    MEMORY_THRESHOLD = 1 * 1024**3  # 1GB
    if data_array.nbytes > MEMORY_THRESHOLD:
        log_message("Using memory-efficient data generator due to large dataset size")
        train_data = DataGenerator(
            data_array[train_idx], batch_size=batch_size, shuffle=True
        )
        val_data = DataGenerator(
            data_array[val_idx], batch_size=batch_size, shuffle=False
        )
        steps_per_epoch = len(train_data)
        validation_steps = len(val_data)
        validation_data = val_data
    else:
        train_data = data_array[train_idx]
        val_data = data_array[val_idx]
        steps_per_epoch = None
        validation_steps = None
        validation_data = (val_data, val_data)

    # -------------------------
    # 2. Basic input validation
    # -------------------------

    input_dim = data_array.shape[1]  # Use numpy array shape

    if n_components >= input_dim:
        log_message(
            f"The input FASTA file has too few sequences to use autoencoder for the dimensionality reduction of kmer counts. n_components ({n_components}) must be less than input dimension ({input_dim}).",
            level="warning",
        )
        return None

    check_model_architecture(input_dim, n_components)

    if nr_of_epochs < 1:
        sys.stderr.write(
            "Warning: Invalid number of epochs. Setting to default of 100.\n"
        )
        nr_of_epochs = 100

    valid_activations = ["relu", "tanh", "sigmoid", "linear", "selu"]
    if activation_mode not in valid_activations:
        sys.stderr.write(
            f"Warning: Invalid activation mode '{activation_mode}'. Using 'relu'\n"
        )
        activation_mode = "relu"

    # -------------------------
    # 3. Build the VAE network
    # -------------------------

    # Encoder with regularization

    try:
        # Create VAE model with Weber's architecture
        vae = VAE(input_dim=input_dim, n_components=n_components)

        # Compile model
        vae.compile(
            optimizer=tf.keras.optimizers.Adam(
                learning_rate=learning_rate, clipnorm=1.0
            )  # Add gradient norm clipping
        )

    except Exception as e:
        log_message(f"Failed to create VAE model: {str(e)}", level="error")
        return None, None if return_umap else None

    # -------------------------
    # 4. Progress callback
    # -------------------------

    class ProgressCallback(tf.keras.callbacks.Callback):
        def __init__(self, total_epochs):
            super().__init__()
            self.total_epochs = total_epochs

        def on_train_begin(self, logs=None):
            log_message(
                f"Starting autoencoder training for {self.total_epochs} epochs",
                level="info",
            )

        def on_epoch_end(self, epoch, logs=None):
            if epoch % 5 == 0:  # Report every 5 epochs
                log_message(
                    f"Epoch {epoch + 1}/{self.total_epochs}: "
                    f"loss = {logs['loss']:.4f}, "
                    f"val_loss = {logs.get('val_loss', 0):.4f}, "
                    f"reconstruction_loss = {logs.get('reconstruction_loss', 0):.4f}, "
                    f"kl_loss = {logs.get('kl_loss', 0):.4f}",
                    level="info",
                )

        def on_train_end(self, logs=None):
            log_message(
                f"Training completed. Final loss: {logs['loss']:.4f}", level="info"
            )

    # -------------------------
    # 5. Warm-up callback
    # -------------------------

    class WarmUpCosineDecay(tf.keras.callbacks.Callback):
        def __init__(self, total_epochs, warmup_epochs=10, base_lr=0.001, min_lr=1e-5):
            super().__init__()
            self.total_epochs = total_epochs
            self.warmup_epochs = warmup_epochs
            self.base_lr = base_lr
            self.min_lr = min_lr
            self.cycle_length = 15  # Increased from 10

        def on_epoch_begin(self, epoch, logs=None):
            if epoch < self.warmup_epochs:
                # Smoother warmup using quadratic growth
                progress = epoch / self.warmup_epochs
                lr = self.base_lr * (progress**2)
            else:
                # Gentler decay
                cycle = (epoch - self.warmup_epochs) // self.cycle_length
                cycle_epoch = (epoch - self.warmup_epochs) % self.cycle_length
                progress = cycle_epoch / self.cycle_length
                lr = self.min_lr + 0.5 * (self.base_lr - self.min_lr) * (
                    1 + np.cos(np.pi * progress)
                ) * (
                    0.9**cycle
                )  # Changed from 0.8

            tf.keras.backend.set_value(self.model.optimizer.lr, lr)

    # -------------------------
    # 6. Dynamic KL weight adjustment callback
    # -------------------------

    class BetaScheduler(tf.keras.callbacks.Callback):
        def __init__(
            self, beta_final=0.0005, ramp_start=20, ramp_length=40
        ):  # Modified values
            super().__init__()
            self.beta_final = beta_final
            self.ramp_start = ramp_start
            self.ramp_length = ramp_length

        def on_epoch_begin(self, epoch, logs=None):
            if epoch < self.ramp_start:
                beta = 0.00001  # Start with even smaller value
            elif epoch < (self.ramp_start + self.ramp_length):
                # More gradual sigmoid ramp
                progress = (epoch - self.ramp_start) / self.ramp_length
                beta = self.beta_final / (
                    1 + np.exp(-5 * (progress - 0.5))
                )  # Reduced slope
            else:
                beta = self.beta_final
            self.model.beta.assign(beta)

    # -------------------------
    # 7. Model performance monitoring
    # -------------------------

    class ModelPerformanceMonitor(tf.keras.callbacks.Callback):
        """
        Monitor various model performance metrics during training.

        Parameters:
            check_frequency (int): How often to check metrics (in epochs)

        Monitors:
            - Reconstruction loss stagnation
            - KL divergence collapse
            - Training stability
            - Gradient norms
        """

        def __init__(self, check_frequency=5):
            super().__init__()
            self.check_frequency = check_frequency
            self.reconstruction_losses = []
            self.kl_losses = []
            self.grad_norms = []

        def on_epoch_end(self, epoch, logs=None):
            if epoch % self.check_frequency != 0:
                return

            # Store losses
            self.reconstruction_losses.append(logs.get("reconstruction_loss", 0))
            self.kl_losses.append(logs.get("kl_loss", 0))

            # Monitor training metrics
            if len(self.reconstruction_losses) > 1:
                # Check for stuck reconstruction
                rec_diff = abs(
                    self.reconstruction_losses[-1] - self.reconstruction_losses[-2]
                )
                if rec_diff < 1e-7:
                    log_message(
                        "Warning: Reconstruction loss has stagnated "
                        f"(change: {rec_diff:.2e}). "
                        "Consider adjusting learning rate or model capacity.",
                        level="warning",
                    )

                # Check for KL collapse
                if self.kl_losses[-1] < 1e-4:
                    log_message(
                        f"Warning: Very low KL divergence detected ({self.kl_losses[-1]:.2e}). "
                        "Model might be experiencing posterior collapse.",
                        level="warning",
                    )

                # Check for unstable KL divergence
                if self.kl_losses[-1] > 10 * self.kl_losses[-2]:
                    log_message(
                        f"Warning: Unstable KL divergence detected "
                        f"(ratio: {self.kl_losses[-1] / self.kl_losses[-2]:.2f}). "
                        "Consider reducing learning rate.",
                        level="warning",
                    )

                # Check for vanishing gradients
                try:
                    with tf.GradientTape() as tape:
                        # Get a batch of data
                        if hasattr(self.model, "train_data"):
                            batch = next(iter(self.model.train_data))
                        else:
                            # If no specific training data attached, use a small random batch
                            batch = tf.random.normal([32, self.model.input_shape[1]])

                        # Compute loss
                        total_loss, _, _ = self.model.compute_loss(batch)

                    # Get gradients
                    gradients = tape.gradient(
                        total_loss, self.model.trainable_variables
                    )
                    grad_norm = tf.linalg.global_norm(gradients)
                    self.grad_norms.append(grad_norm.numpy())

                    # Check gradient norm
                    if grad_norm < 1e-7:
                        log_message(
                            f"Warning: Very small gradient norm detected ({grad_norm:.2e}). "
                            "Model might be experiencing vanishing gradients.",
                            level="warning",
                        )
                    elif grad_norm > 10.0:
                        log_message(
                            f"Warning: Large gradient norm detected ({grad_norm:.2e}). "
                            "Consider reducing learning rate or increasing gradient clipping threshold.",
                            level="warning",
                        )

                except Exception as e:
                    log_message(
                        f"Note: Could not compute gradient metrics: {str(e)}",
                        level="info",
                    )

            # Log current state
            log_message(
                f"Epoch {epoch + 1} metrics - "
                f"Reconstruction loss: {self.reconstruction_losses[-1]:.4f}, "
                f"KL loss: {self.kl_losses[-1]:.4f}",
                level="info",
            )

    # -------------------------
    # 8. Compile + training
    # -------------------------

    learning_rate = 0.001  # Base learning rate

    if nr_of_epochs > 100:
        log_message(
            f"Training with {nr_of_epochs} epochs may take a long time. "
            "Consider reducing if training time is a concern. "
            "Early stopping will prevent overtraining.",
            level="warning",
        )

    # Update the callbacks list:
    callbacks = [
        WarmUpCosineDecay(
            total_epochs=nr_of_epochs,
            warmup_epochs=10,
            base_lr=learning_rate,
            min_lr=1e-5,
        ),
        BetaScheduler(
            beta_final=0.0005,  # Match class default
            ramp_start=20,  # Match class default
            ramp_length=40,  # Match class default
        ),
        ProgressCallback(nr_of_epochs),
        tf.keras.callbacks.EarlyStopping(
            monitor="val_loss",
            patience=30,
            restore_best_weights=True,
            verbose=1,
            min_delta=1e-4,  # Increased patience
        ),
        tf.keras.callbacks.ReduceLROnPlateau(
            monitor="val_loss",
            factor=0.5,  # More gradual reduction
            patience=10,  # Increased patience
            min_lr=1e-6,
            verbose=1,
            mode="min",
        ),
        tf.keras.callbacks.LambdaCallback(
            on_epoch_end=lambda epoch, logs: (
                log_message(
                    f"Epoch {epoch + 1}: learning rate = {K.get_value(vae.optimizer.learning_rate):.6f}",
                    level="info",
                )
                if epoch % 5 == 0
                else None
            )
        ),
        ModelPerformanceMonitor(check_frequency=5),
        tf.keras.callbacks.CSVLogger(
            os.path.join(out_folder, "training_log.csv"), separator=",", append=False
        ),
    ]

    # Modified training data handling
    if steps_per_epoch is None:
        # For small datasets, use simple numpy arrays
        validation_data = val_data  # Just the data, not a tuple
    else:
        # For large datasets, use generator
        validation_data = val_data

    # Memory cleanup
    gc.collect()

    # Training call
    history = vae.fit(
        train_data,
        validation_data=validation_data if validation_steps else (val_data, val_data),
        epochs=nr_of_epochs,
        batch_size=batch_size if steps_per_epoch is None else None,
        callbacks=callbacks,
        shuffle=True,
        steps_per_epoch=steps_per_epoch,
        validation_steps=validation_steps,
        verbose=0,
    )

    # Memory cleanup after training
    gc.collect()

    # -------------------------
    # 9. Save training plots
    # -------------------------
    if out_folder is not None:
        monitor_vae_training(history, out_folder)

    # -------------------------
    # 10. Build encoder for embeddings
    # -------------------------
    encoder = vae.encoder

    # For prediction, we might still want to use batches for large datasets
    embedding = None
    if data_array.shape[0] > MEMORY_THRESHOLD:
        embeddings = []
        for i in range(0, len(data_array), batch_size):
            batch = data_array[i : i + batch_size]
            # Get only z_mean from the encoder output
            z_mean, _ = encoder.predict(batch)
            embeddings.append(z_mean)
        embedding = np.vstack(embeddings)
    else:
        # Get only z_mean from the encoder output
        z_mean, _ = encoder.predict(data_array)
        embedding = z_mean

    if np.any(np.isnan(embedding)):
        n_nans = np.isnan(embedding).sum()
        log_message(
            f"Found {n_nans} NaN values in embedding. "
            "These will be replaced with zeros. This might indicate "
            "training instability.",
            level="warning",
        )
        embedding = np.nan_to_num(embedding, nan=0.0)

    if np.any(np.isnan(embedding)):
        log_message("Unable to fix NaN values in embedding", level="error")
        return None, None if return_umap else None

    # -------------------------
    # 11. Running UMAP with autoencoder's embedding
    # -------------------------

    umap_embedding = None
    if return_umap:

        # UMAP visualisation of the latent space
        if optimise_n_neighbors:
            log_message(
                "Optimising UMAP parameters for latent space visualisation...",
                level="info",
            )
            best_n_neighbors, opt_results = find_best_n_neighbors_setting(
                embedding,
                method_name="umap",
                n_neighbors_range=[5, 10, 15, 20, 30],
                min_cluster_sizes=[5, 10, 15],
            )

            # Save optimisation results
            if out_folder:
                results_path = os.path.join(
                    out_folder, "autoencoder_umap_optimisation.txt"
                )
                with open(results_path, "w") as f:
                    f.write("n_neighbors,score\n")
                    for n, score in opt_results.items():
                        f.write(f"{n},{score:.4f}\n")

            log_message(
                f"Best n_neighbors for UMAP visualisation: {best_n_neighbors}",
                level="info",
            )
            umap_embedding = umap.UMAP(
                n_neighbors=best_n_neighbors, n_components=5
            ).fit_transform(embedding)
        else:
            # Use default UMAP parameters
            umap_embedding = umap.UMAP(n_components=5).fit_transform(embedding)

        # Save UMAP visualisation
        if out_folder is not None:
            save_embedding_results(
                embedding=umap_embedding,
                original_data=data_array,  # Still use original data for comparison
                out_folder=out_folder,
                method_name="UMAP_of_VAE",
                seq_names=seq_names,
            )

    if return_umap:
        return embedding, umap_embedding
    return embedding


def monitor_vae_training(history, out_folder):
    """Enhanced training visualisation with additional metrics."""
    plt.figure(figsize=(15, 10))

    # Plot total loss
    plt.subplot(2, 2, 1)
    plt.plot(history.history["loss"], label="Train")
    plt.plot(history.history["val_loss"], label="Validation")
    plt.title("Total Loss")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.legend()

    # Plot reconstruction loss
    plt.subplot(2, 2, 2)
    plt.plot(history.history["reconstruction_loss"], label="Train")
    if "val_reconstruction_loss" in history.history:
        plt.plot(history.history["val_reconstruction_loss"], label="Validation")
    plt.title("Reconstruction Loss")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.legend()

    # Plot KL loss
    plt.subplot(2, 2, 3)
    plt.plot(history.history["kl_loss"], label="Train")
    if "val_kl_loss" in history.history:
        plt.plot(history.history["val_kl_loss"], label="Validation")
    plt.title("KL Divergence Loss")
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.legend()

    # Plot learning rate
    plt.subplot(2, 2, 4)
    if "lr" in history.history:
        plt.plot(history.history["lr"], label="Learning Rate")
        plt.title("Learning Rate over Time")
        plt.xlabel("Epoch")
        plt.ylabel("Learning Rate")
        plt.yscale("log")

    plt.tight_layout()
    plt.savefig(os.path.join(out_folder, "vae_training_history.png"))
    plt.close()


def visualise_embedding(
    embedding, metric_values=None, out_folder=None, title="Embedding"
):
    """
    Create comprehensive visualisation of embedding results with robust error handling.
    """
    if embedding.shape[1] < 2:
        log_message("Cannot visualise 1D embedding", level="warning")
        return

    plt.figure(figsize=(15, 10))

    # Main scatter plot
    plt.subplot(2, 2, 1)
    if metric_values is not None and len(metric_values) == embedding.shape[0]:
        scatter = plt.scatter(
            embedding[:, 0], embedding[:, 1], c=metric_values, cmap="viridis", alpha=0.6
        )
        plt.colorbar(scatter, label="Quality score")
    else:
        plt.scatter(embedding[:, 0], embedding[:, 1], alpha=0.6)

    plt.title(f"{title} - Main View")
    plt.xlabel("Component 1")
    plt.ylabel("Component 2")

    # 2D histogram (density) plot
    plt.subplot(2, 2, 2)
    try:
        hist2d = plt.hist2d(embedding[:, 0], embedding[:, 1], bins=50, cmap="viridis")
        plt.colorbar(hist2d[3], label="Count")
    except Exception as e:
        log_message(f"Could not create 2D histogram: {str(e)}", level="warning")
        plt.text(0.5, 0.5, "Histogram unavailable", ha="center", va="center")

    plt.title(f"{title} - Density View")
    plt.xlabel("Component 1")
    plt.ylabel("Component 2")

    # Distribution plots with robust KDE
    for idx, pos in enumerate([(2, 2, 3), (2, 2, 4)]):
        plt.subplot(*pos)
        try:
            # Simple histogram
            plt.hist(embedding[:, idx], bins=50, density=True, alpha=0.7)

            # Robust KDE estimation
            try:
                from scipy.stats import gaussian_kde

                # Add small noise to prevent singular matrix
                data = embedding[:, idx] + np.random.normal(
                    0, 1e-10, embedding[:, idx].shape
                )
                kde = gaussian_kde(data, bw_method="silverman")
                x_range = np.linspace(data.min(), data.max(), 200)
                plt.plot(x_range, kde(x_range))
            except Exception as e:
                log_message(
                    f"KDE failed for component {idx+1}: {str(e)}", level="warning"
                )
                # Continue without KDE line
        except Exception as e:
            log_message(
                f"Could not create distribution plot for component {idx+1}: {str(e)}",
                level="warning",
            )
            plt.text(0.5, 0.5, "Distribution unavailable", ha="center", va="center")

        plt.title(f"Distribution of Component {idx+1}")
        plt.xlabel("Value")
        plt.ylabel("Density")

    plt.tight_layout()

    if out_folder:
        try:
            plt.savefig(
                os.path.join(
                    out_folder, f"{title.lower().replace(' ', '_')}_visualisation.png"
                ),
                bbox_inches="tight",
                dpi=300,
            )
        except Exception as e:
            log_message(f"Failed to save visualisation: {str(e)}", level="error")
        finally:
            plt.close()


def calculate_embedding_quality(original_data, embedding, n_neighbors=15):
    """
    Calculate trustworthiness and continuity of an embedding.

    Parameters:
        original_data: array, shape (n_samples, n_features)
            Original high-dimensional data
        embedding: array, shape (n_samples, n_components)
            Low-dimensional embedding
        n_neighbors: int
            Number of neighbors to consider

    Returns:
        dict with trustworthiness and continuity scores
    """
    if (
        original_data.shape[0] < 3
    ):  # Absolute minimum for meaningful neighborhood analysis
        log_message("Dataset too small for meaningful quality metrics", level="warning")
        return {"trustworthiness": float("nan"), "continuity": float("nan")}

    n_samples = original_data.shape[0]

    # Ensure n_neighbors is valid
    n_neighbors = min(n_neighbors, n_samples - 1)

    # Fit nearest neighbors on both spaces
    nn_orig = NearestNeighbors(n_neighbors=n_neighbors + 1)  # +1 for self
    nn_emb = NearestNeighbors(n_neighbors=n_neighbors + 1)

    nn_orig.fit(original_data)
    nn_emb.fit(embedding)

    # Get nearest neighbors (indices only)
    orig_ind = nn_orig.kneighbors(return_distance=False)
    emb_ind = nn_emb.kneighbors(return_distance=False)

    # Remove self-references (first column)
    orig_ind = orig_ind[:, 1:]
    emb_ind = emb_ind[:, 1:]

    # Calculate trustworthiness
    trustworthiness = 0
    for i in range(n_samples):
        # Find points that are in k nearest neighbors in embedding but not in original
        emb_neighs = set(emb_ind[i])
        orig_neighs = set(orig_ind[i])
        new_neighs = emb_neighs - orig_neighs

        if new_neighs:
            # Find rank in original space
            for j in new_neighs:
                rank = (
                    np.where(orig_ind[i] == j)[0][0] + 1
                    if j in orig_ind[i]
                    else n_neighbors + 1
                )
                trustworthiness += rank - n_neighbors

    # Calculate continuity
    continuity = 0
    for i in range(n_samples):
        # Find points that are in k nearest neighbors in original but not in embedding
        orig_neighs = set(orig_ind[i])
        emb_neighs = set(emb_ind[i])
        lost_neighs = orig_neighs - emb_neighs

        if lost_neighs:
            # Find rank in embedding space
            for j in lost_neighs:
                rank = (
                    np.where(emb_ind[i] == j)[0][0] + 1
                    if j in emb_ind[i]
                    else n_neighbors + 1
                )
                continuity += rank - n_neighbors

    # Normalise
    n = n_samples * n_neighbors * (2 * n_samples - 3 * n_neighbors - 1) / 2
    trustworthiness = 1 - trustworthiness / n
    continuity = 1 - continuity / n

    return {"trustworthiness": float(trustworthiness), "continuity": float(continuity)}


def save_embedding_results(
    embedding, original_data, out_folder, method_name, seq_names=None
):
    """
    Save embedding results including visualisations and quality metrics.
    """
    if not os.path.exists(out_folder):
        try:
            os.makedirs(out_folder)
        except OSError as e:
            log_message(f"Error creating output directory: {e}", level="error")
            return

    # Check embedding quality
    if np.any(np.isnan(embedding)):
        log_message(
            "Embedding contains NaN values - saving will be limited", level="warning"
        )
        return

    if np.all(embedding == 0):
        log_message(
            "Embedding contains all zeros - visualisation may be uninformative",
            level="warning",
        )

    # Calculate quality metrics
    try:
        quality_metrics = calculate_embedding_quality(original_data, embedding)

        # Save metrics
        with open(os.path.join(out_folder, f"{method_name}_metrics.txt"), "w") as f:
            for metric, value in quality_metrics.items():
                f.write(f"{metric}: {value:.4f}\n")
    except Exception as e:
        log_message(f"Failed to calculate quality metrics: {str(e)}", level="warning")

    # Create visualisation with error handling
    try:
        visualise_embedding(
            embedding, metric_values=None, out_folder=out_folder, title=method_name
        )
    except Exception as e:
        log_message(f"Failed to create visualisation: {str(e)}", level="error")

    # Save embedding coordinates
    if seq_names is not None and "_raw" in method_name:
        try:
            embedding_df = embedding_to_dataframe(embedding, seq_names, method_name)
            embedding_df.to_csv(
                os.path.join(out_folder, f"{method_name}_coordinates.csv"), index=False
            )
        except Exception as e:
            log_message(
                f"Failed to save embedding coordinates: {str(e)}", level="error"
            )

    # Save clustering metrics
    if seq_names is not None:
        # Extract organism labels from seq_names
        labels = [name.split("_")[0] for name in seq_names]

        # Calculate clustering metrics
        metrics = evaluate_clustering(embedding, labels)

        # Save metrics
        with open(
            os.path.join(out_folder, f"{method_name}_clustering_metrics.txt"), "w"
        ) as f:
            for metric, value in metrics.items():
                f.write(f"{metric}: {value:.4f}\n")


def run_tsne(df, n_components, n_neighbors):
    """
    Function for running t-SNE for dimensionality reduction of kmer counts
    """
    print_timestamp("t-SNE")
    df_row_count = df.shape[0]
    if df_row_count < 2:
        log_message(
            "Cannot run t-SNE because the kmer counts dataframe has only {} rows\n".format(
                df_row_count
            ),
            level="warning",
        )
        return None

    if n_components > 3:
        # log_message(f"Warning: t-SNE typically performs best with 2 or 3 components. Results with {n_components} components may be less reliable.", level="warning")
        log_message(
            f"The user-defined n_components value {n_components} does not work with the barnes_hut algorithm that t-SNE uses, as it relies on quad-tree or oct-tree. Reducing n_components to 3 for t-SNE, in order to be able to run it",
            level="warning",
        )
        n_components = 3

    if n_components <= 0:
        log_message("n_components must be positive", level="error")
        return None

    # Use n_neighbors_setting as perplexity, but ensure it's valid
    perplexity = min(
        n_neighbors, df_row_count / 3
    )  # t-SNE typically expects perplexity < N/3
    if perplexity < 1:
        perplexity = 1
    log_message(f"Using t-SNE perplexity value of {perplexity}", level="warning")

    # Add check for minimum number of samples needed
    min_samples_needed = n_components + 1
    if df_row_count < min_samples_needed:
        log_message(
            f"Cannot run t-SNE with {n_components} components because the dataset has only {df_row_count} samples. Need at least {min_samples_needed} samples.",
            level="warning",
        )
        return None

    embedding = manifold.TSNE(
        n_components=n_components, perplexity=perplexity
    ).fit_transform(df)
    return embedding


def get_min_neighbors_for_lle(method, n_components):
    """
    Calculate minimum required neighbors for different LLE methods
    """
    try:
        if method == "hessian":
            return max(
                int((n_components * (n_components + 3) / 2) + 1), n_components + 1
            )
        elif method == "standard":
            return max(n_components + 1, 2)  # At least 2 neighbors
        elif method == "modified":
            return max(n_components + 2, 2)  # At least 2 neighbors
        else:
            return max(n_components + 1, 2)  # Default fallback
    except Exception as e:
        log_message(f"Error calculating minimum neighbors: {str(e)}", level="warning")
        return n_components + 2  # Safe fallback


def run_lle(df_preprocessed, method, n_neighbors_setting, n_components):
    """
    Run LLE with proper error handling for all variants
    """
    method_names = {
        "standard": "LLE standard",
        "hessian": "LLE hessian",
        "modified": "LLE modified",
    }

    #
    if n_neighbors_setting is None:
        n_neighbors_setting = n_components + 2  # Safe default
        log_message(
            f"No valid n_neighbors setting provided for {method} LLE. "
            f"Using default value: {n_neighbors_setting}",
            level="warning",
        )

    embedding_title = method_names.get(method, f"LLE {method}")
    print_timestamp(embedding_title)

    # Calculate minimum required neighbors
    min_neighbors = get_min_neighbors_for_lle(method, n_components)

    if n_neighbors_setting < min_neighbors:
        n_neighbors_setting = min_neighbors
        log_message(
            f"Adjusted n_neighbors to minimum required value ({min_neighbors}) "
            f"for {embedding_title}",
            level="warning",
        )

    try:
        # Try with dense solver first
        embedding = manifold.LocallyLinearEmbedding(
            n_neighbors=n_neighbors_setting,
            n_components=n_components,
            eigen_solver="dense",
            method=method,
        ).fit_transform(df_preprocessed)
    except Exception as e:
        log_message(
            f"Error running {embedding_title} with dense solver: {str(e)}. "
            "Trying with arpack solver...",
            level="warning",
        )
        try:
            # Try with arpack solver as fallback
            embedding = manifold.LocallyLinearEmbedding(
                n_neighbors=n_neighbors_setting,
                n_components=n_components,
                eigen_solver="arpack",
                method=method,
            ).fit_transform(df_preprocessed)
        except Exception as e:
            log_message(
                f"Error running {embedding_title} with arpack solver: {str(e)}. "
                "This method might not be suitable for your data.",
                level="error",
            )
            embedding = None

    return embedding, embedding_title


def check_minimum_samples(n_samples, n_features, method, n_components):
    """
    Check if there are enough samples for the selected method.
    Returns (bool, str, int): (is_valid, message, n_components)
    """

    method_requirements = {
        "pca": {"min_samples": 2, "message": "PCA requires at least 2 samples"},
        "random_trees": {
            "min_samples": 2,
            "message": "Random Trees embedding requires at least 2 samples",
        },
        "umap": {
            "min_samples": 10,
            "message": "UMAP works best with at least 10 samples",
        },
        "t-sne": {"min_samples": 5, "message": "t-SNE requires at least 5 samples"},
        "isomap": {"min_samples": 5, "message": "Isomap requires at least 5 samples"},
        "lle_standard": {
            "min_samples": max(n_components + 2, 5),
            "message": f"Standard LLE requires at least {max(n_components + 2, 5)} samples",
        },
        "lle_hessian": {
            "min_samples": max(n_components * (n_components + 3) // 2 + 1, 5),
            "message": f"Hessian LLE requires at least {max(n_components * (n_components + 3) // 2 + 1, 5)} samples",
        },
        "lle_modified": {
            "min_samples": max(n_components + 2, 5),
            "message": f"Modified LLE requires at least {max(n_components + 2, 5)} samples",
        },
        "mds": {"min_samples": 4, "message": "MDS requires at least 4 samples"},
        "se": {
            "min_samples": 4,
            "message": "Spectral Embedding requires at least 4 samples",
        },
        "kernel_pca": {
            "min_samples": 2,
            "message": "Kernel PCA requires at least 2 samples",
        },
        "pca_svd": {
            "min_samples": 2,
            "message": "PCA with SVD requires at least 2 samples",
        },
        "nmf": {"min_samples": 2, "message": "NMF requires at least 2 samples"},
        "autoencoder_relu": {
            "min_samples": 10,
            "message": "Autoencoder requires at least 10 samples for meaningful results",
        },
    }

    # Add common checks for all autoencoder variants
    if method.startswith("autoencoder_"):
        method_requirements[method] = {
            "min_samples": 10,
            "message": "Autoencoder requires at least 10 samples for meaningful results",
        }

    if method not in method_requirements:
        return True, ""  # Unknown method, assume it's valid

    req = method_requirements[method]
    is_valid = n_samples >= req["min_samples"]

    # Additional check for n_components
    if is_valid and n_components >= n_samples:
        is_valid = False
        req["message"] = (
            f"Number of components ({n_components}) must be less than number of samples ({n_samples})"
        )

    # Check for feature count
    if is_valid and n_features < 2:
        is_valid = False
        req["message"] = "At least 2 features are required"

    return is_valid, req["message"], n_components


def run_dim_reduction(
    df,
    seq_lengths,
    seq_names,
    selected_method,
    n_neighbors_setting=-1,
    autoencoder_epochs_count=-1,
    out_folder=None,
    n_components=5,
    optimise_n_neighbors=False,
):

    embedding = None
    embedding_title = None

    n_samples, n_features = df.shape

    # Adjust n_components before anything else
    if n_components >= n_samples:
        n_components = max(1, n_samples - 1)
        log_message(
            f"Adjusting n_components from {n_components + 1} to {n_components} "
            f"for {selected_method} due to small sample size",
            level="warning",
        )

    # Check if we have enough samples for the selected method
    is_valid, message, _ = check_minimum_samples(
        n_samples, n_features, selected_method, n_components
    )

    if not is_valid:
        log_message(f"Skipping {selected_method}: {message}", level="warning")
        return None, None, False

    methods_with_neighbors = {
        "umap",
        "isomap",
        "lle_standard",
        "lle_hessian",
        "lle_modified",
        "se",
    }

    # Validate n_components at the start
    n_samples, n_features = df.shape

    if n_components <= 0:
        log_message(f"Invalid n_components value ({n_components})", level="error")
        return None, None

    if optimise_n_neighbors and selected_method in methods_with_neighbors:
        log_message(f"Optimising parameters for {selected_method}...", level="info")
        try:
            best_n_neighbors, opt_results = find_best_n_neighbors_setting(
                df.values, selected_method
            )

            # Save optimisation results
            if out_folder:
                results_path = os.path.join(
                    out_folder, f"{selected_method}_n_neighbors_optimisation.txt"
                )
                with open(results_path, "w") as f:
                    f.write("n_neighbors,score\n")
                    for n, score in opt_results.items():
                        f.write(f"{n},{score:.4f}\n")

            if best_n_neighbors is not None:
                n_neighbors_setting = best_n_neighbors
                log_message(
                    f"Best n_neighbors for {selected_method}: {best_n_neighbors}",
                    level="info",
                )
            else:
                # Use default or provided n_neighbors if optimisation fails
                log_message(
                    f"n_neighbors optimisation failed for {selected_method}. "
                    f"Using provided/default n_neighbors={n_neighbors_setting}",
                    level="warning",
                )
        except Exception as e:
            log_message(
                f"Error during n_neighbors parameter optimisation for {selected_method}: {str(e)}. "
                f"Using provided/default n_neighbors={n_neighbors_setting}",
                level="warning",
            )

    if selected_method in methods_with_neighbors:
        if n_neighbors_setting is None or n_neighbors_setting < 1:
            n_neighbors_setting = n_components + 2  # Safe default
            log_message(
                f"Invalid n_neighbors setting for {selected_method}. "
                f"Using default value: {n_neighbors_setting}",
                level="warning",
            )
        n_neighbors_setting = validate_n_neighbors(
            n_neighbors_setting, n_samples, selected_method
        )

    # Preprocessing based on the method
    if selected_method in [
        "pca",
        "lle_standard",
        "lle_hessian",
        "lle_modified",
        "kernel_pca",
        "pca_svd",
        "se",
    ]:
        # Standard scaling ensures features contribute equally in PCA/LLE
        df_preprocessed = preprocess_data(
            df, seq_lengths, scaling_method="standard", use_tf_idf=False
        )
    elif selected_method in ["umap", "t-sne", "isomap", "mds"]:
        df_preprocessed = preprocess_data(
            df, seq_lengths, scaling_method="log", use_tf_idf=False
        )
    elif selected_method in ["random_trees", "nmf"]:
        # df_preprocessed = df.values  # Pass raw values because these methods don't require scaling
        df_preprocessed = preprocess_data(
            df, seq_lengths, scaling_method="log", use_tf_idf=False
        )
    elif selected_method in [
        "autoencoder_sigmoid",
        "autoencoder_linear",
        "autoencoder_tanh",
        "autoencoder_selu",
        "autoencoder_relu",
    ]:
        df_preprocessed = df.values  # Pass raw values; autoencoder handles scaling
    else:
        raise ValueError(f"Unknown dimensionality reduction method: {selected_method}")

    if selected_method == "pca":
        embedding_title = "PCA"
        print_timestamp(embedding_title)
        pca = decomposition.PCA(n_components=n_components)
        embedding = pca.fit_transform(df_preprocessed)

    elif selected_method == "umap":
        embedding_title = "UMAP"
        print_timestamp(embedding_title)
        reducer = umap.UMAP(
            random_state=123,
            n_neighbors=n_neighbors_setting,
            min_dist=0.1,
            n_components=n_components,
            metric="euclidean",
        )
        embedding = reducer.fit_transform(df_preprocessed)

    elif selected_method == "t-sne":
        embedding_title = "t-SNE"
        embedding = run_tsne(
            df_preprocessed, n_neighbors=n_neighbors_setting, n_components=n_components
        )

    elif selected_method == "isomap":
        embedding_title = "isomap"
        print_timestamp(embedding_title)
        embedding = manifold.Isomap(
            n_neighbors=n_neighbors_setting, n_components=n_components
        ).fit_transform(df_preprocessed)

    elif selected_method == "lle_standard":
        embedding, embedding_title = run_lle(
            df_preprocessed, "standard", n_neighbors_setting, n_components
        )
    elif selected_method == "lle_hessian":
        embedding, embedding_title = run_lle(
            df_preprocessed, "hessian", n_neighbors_setting, n_components
        )
    elif selected_method == "lle_modified":
        embedding, embedding_title = run_lle(
            df_preprocessed, "modified", n_neighbors_setting, n_components
        )

    elif selected_method == "mds":
        embedding_title = "MDS"
        print_timestamp(embedding_title)
        # MDS can't produce more components than n_samples - 1
        max_possible_components = min(df.shape[0] - 1, df.shape[1])
        if n_components > max_possible_components:
            log_message(
                f"Warning: requested {n_components} components but MDS can only compute {max_possible_components} components with this dataset. Reducing n_components to {max_possible_components}",
                level="warning",
            )
            n_components = max_possible_components
        embedding = manifold.MDS(
            max_iter=100, n_init=1, n_components=n_components
        ).fit_transform(df_preprocessed)

    elif selected_method == "se":
        embedding_title = "SE"
        print_timestamp(embedding_title)
        embedding = manifold.SpectralEmbedding(
            n_neighbors=n_neighbors_setting, n_components=n_components
        ).fit_transform(df_preprocessed)

    elif selected_method == "random_trees":
        # https://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html#sphx-glr-auto-examples-manifold-plot-lle-digits-py
        embedding_title = "random_trees"
        print_timestamp(embedding_title)
        hasher = ensemble.RandomTreesEmbedding(
            n_estimators=200, random_state=0, max_depth=5
        )
        x_transformed = hasher.fit_transform(df_preprocessed)
        pca = decomposition.TruncatedSVD(n_components=n_components)
        embedding = pca.fit_transform(x_transformed)

    elif selected_method == "kernel_pca":
        # https://www.tutorialspoint.com/scikit_learn/scikit_learn_dimensionality_reduction_using_pca.htm
        embedding_title = "KernelPCA"
        print_timestamp(embedding_title)
        kernel_pca = decomposition.KernelPCA(
            n_components=n_components, kernel="sigmoid"
        )
        embedding = kernel_pca.fit_transform(df_preprocessed)

    elif selected_method == "pca_svd":
        # https://www.tutorialspoint.com/scikit_learn/scikit_learn_dimensionality_reduction_using_pca.htm
        embedding_title = "PCA with SVD solver"
        print_timestamp(embedding_title)
        pca = decomposition.PCA(n_components=n_components, svd_solver="randomized")
        embedding = pca.fit_transform(df_preprocessed)

    elif selected_method == "nmf":
        # https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html
        embedding_title = "Non-Negative Matrix Factorization"
        print_timestamp(embedding_title)
        nmf_model = decomposition.NMF(
            n_components=n_components, init="random", random_state=0
        )
        embedding = nmf_model.fit_transform(df_preprocessed)

    elif selected_method.startswith("autoencoder"):
        scaling_method = "log"  # Adjust if needed
        activation_mode = selected_method.split("_")[1]
        embedding_title = f"Autoencoder ({activation_mode})"
        print_timestamp(embedding_title)

        # Get both embeddings
        raw_embedding, umap_embedding = run_autoencoder(
            df,
            seq_lengths,
            autoencoder_epochs_count,
            activation_mode,
            out_folder,
            n_components=n_components,
            scaling_method=scaling_method,
            return_umap=True,
            optimise_n_neighbors=optimise_n_neighbors,
        )

        if raw_embedding is not None:
            # Save raw autoencoder embedding visualisation and coordinates
            save_embedding_results(
                embedding=raw_embedding,
                original_data=df.values,
                out_folder=out_folder,
                method_name=f"{selected_method}_raw",
                seq_names=seq_names,
            )

        # Use UMAP embedding as the main output
        embedding = umap_embedding
        embedding_title = f"Autoencoder_{activation_mode}_UMAP"
        return embedding, embedding_title, True

    return embedding, embedding_title, False


def main(
    kmer_counts_file,
    out_folder,
    selected_methods,
    n_neighbors_setting,
    autoencoder_epochs_count,
    n_components,
    skip_n_neighbors_optimisation,
):
    reset_random_seeds()

    # Invert the skip flag to get optimise_n_neighbors
    optimise_n_neighbors = not skip_n_neighbors_optimisation

    if optimise_n_neighbors:
        log_message(
            "n_neighbors optimisation enabled for applicable methods (umap, isomap, lle_standard, lle_hessian, lle_modified, se)",
            level="info",
        )
    else:
        log_message(
            f"n_neighbors optimisation disabled, using provided/default value: {n_neighbors_setting}",
            level="info",
        )

    df, seq_names, seq_lengths = load_data(kmer_counts_file)

    # Remove constant features (k-mers)
    df, n_constant_removed = remove_constant_features(df)

    # Log the new dimensionality
    if n_constant_removed > 0:
        log_message(
            f"Dataset dimensions after removing constant k-mers: {df.shape[0]} samples x {df.shape[1]} k-mers",
            level="info",
        )

    # Early checks for minimal dataset
    df_row_count = df.shape[0]
    if df_row_count == 0:
        log_message("Empty input file", level="error")
        sys.exit(1)

    if df_row_count == 1:
        log_message(
            "Skipping the dimensionality reduction of kmer counts, as the kmer counts table has only one row",
            level="warning",
        )
        # Generate an empty file to satisfy nextflow expecting a file from script finishing with no file with small output
        with open(
            os.path.join(
                os.getcwd(),
                f"EMPTY_{selected_methods}_kmers_dim_reduction_embeddings.csv",
            ),
            "w",
        ) as empty_file:
            empty_file.write("The kmer counts file is too small for analysis")
        sys.exit(0)

    selected_methods = [m.strip() for m in selected_methods.split(",")]

    if df_row_count < n_components:
        log_message(
            f"Warning: Input data has {df_row_count} sequences. Because of this, the user-specified value of n_components ({n_components}) cannot be used. Changing the value of {n_components} to {df_row_count} to make it fit the low number of sequences in the input FASTA file",
            level="warning",
        )
        n_components = df_row_count

    Path(out_folder).mkdir(parents=True, exist_ok=True)
    check_memory_requirements(df, selected_methods)

    if autoencoder_epochs_count == -1:
        autoencoder_epochs_count = df_row_count

    embeddings_list = list()

    for selected_method in selected_methods:
        embedding, embedding_title, skip_visualisation = run_dim_reduction(
            df,
            seq_lengths,
            seq_names,
            selected_method,
            n_neighbors_setting=n_neighbors_setting,
            autoencoder_epochs_count=autoencoder_epochs_count,
            out_folder=out_folder,
            n_components=n_components,
            optimise_n_neighbors=optimise_n_neighbors,
        )

        if embedding is not None:
            # Save additional visualisation and quality metrics only if not skipping
            if not skip_visualisation:
                save_embedding_results(
                    embedding=embedding,
                    original_data=df,
                    out_folder=out_folder,
                    method_name=selected_method,
                    seq_names=seq_names,
                )

            # Keep original embedding export for downstream compatibility
            embedding_df = embedding_to_dataframe(embedding, seq_names, embedding_title)
            embeddings_list.append(embedding_df)

    if not embeddings_list:
        log_message(
            "No valid embeddings were generated. Check the warnings above.",
            level="error",
        )
        with open(
            f"EMPTY_{selected_methods}_kmers_dim_reduction_embeddings.csv", "w"
        ) as empty_file:
            empty_file.write("NO VALID EMBEDDINGS GENERATED")
        sys.exit(1)

    out_df = reduce(
        lambda left, right: pd.merge(
            left, right, on=["scaff"], how="outer"
        ),  # Merge DataFrames in list
        embeddings_list,
    )
    out_path = os.path.join(out_folder, "kmers_dim_reduction_embeddings.csv")
    out_df.to_csv(out_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-v", "--version", action="version", version="1.0")
    parser.add_argument(
        "kmer_counts_file", type=str, help="Path to input CSV file with kmer counts"
    )
    parser.add_argument(
        "out_folder", type=str, help="Path to folder where output files will be written"
    )
    parser.add_argument(
        "--selected_methods",
        type=str,
        help="Comma separated string with the selected dimensionality reduction methods. It should consist of items from this list: pca,umap,t-sne,isomap,lle_standard,lle_hessian,lle_modified,mds,se,random_trees,kernel_pca,pca_svd,nmf,autoencoder_sigmoid,autoencoder_linear,autoencoder_tanh,autoencoder_selu,autoencoder_relu. Default: PCA",
        default="pca",
    )
    parser.add_argument(
        "--n_neighbors_setting",
        type=int,
        help="n_neighbors parameter default value for the methods that have this parameter (umap, isomap, lle_standard, lle_hessian, lle_modified, se). Default: 13. This is used if n_neighbors optimisation is not used",
        default=13,
    )
    parser.add_argument(
        "--autoencoder_epochs_count",
        type=int,
        help="Autoencoder epochs count (default: assign automatically as 3x number of sequences in input)",
        default=-1,
    )
    parser.add_argument(
        "--n_components",
        type=int,
        default=5,
        help="Number of dimensions to retain in the reduced representation (default: 5)",
    )
    parser.add_argument(
        "--skip_n_neighbors_optimisation",
        action="store_true",
        help="Skip the optimisation of n_neighbors for applicable methods (umap, isomap, lle_standard, lle_hessian, lle_modified, se)",
    )
    args = parser.parse_args()
    main(
        args.kmer_counts_file,
        args.out_folder,
        args.selected_methods,
        args.n_neighbors_setting,
        args.autoencoder_epochs_count,
        args.n_components,
        args.skip_n_neighbors_optimisation,
    )
