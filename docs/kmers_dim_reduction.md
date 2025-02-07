# K-mer Dimensionality Reduction for Taxonomic Analysis

## Introduction

The k-mer dimensionality reduction tool is a component of the Assembly Screen for Cobionts and Contaminants (ASCC) pipeline, designed to enhance taxonomic analysis of genome assemblies. While traditional methods like [BlobToolKit](https://blobtoolkit.genomehubs.org/) visualise sequences using GC content and coverage metrics (see [example GC-coverage plot](https://blobtoolkit.genomehubs.org/wp-content/uploads/2019/07/ACVV01.blob_.square-1024x1024.png)), this tool provides complementary visualisations based on sequence composition through k-mer analysis.

### Purpose and Applications

The tool creates dimensionality reduction embeddings of k-mer counts, enabling the visualisation of taxonomic relationships between sequences in an assembly. The fundamental principle is that sequences sharing taxonomic origin tend to cluster together in these embeddings, with organellar sequences typically appearing distinct from their associated chromosomal sequences.

This approach is particularly valuable when:

- Analysing assemblies with complex species composition
- Working with environmental samples
- Examining assemblies containing 3-10 distinct taxonomic groups
- Processing datasets with 100 to several thousand sequences

### Integration with [BlobToolKit](https://blobtoolkit.genomehubs.org/)

The tool seamlessly integrates with BlobToolKit, allowing users to:

- Add custom dimensionality reduction coordinates as variables in BlobToolKit datasets
- Interactively explore and select sequence clusters
- Generate lists of sequences belonging to specific clusters
- Complement traditional GC-coverage plots with alternative taxonomic visualisations

### Data Normalization

Before applying dimensionality reduction methods, the k-mer count data is normalized to a range of 0 to 1. This normalization step:

- Ensures all features contribute equally to the analysis
- Improves numerical stability during computation
- Makes the results more comparable across different sequences
- Helps prevent any single feature from dominating the analysis

## Available Methods

The tool implements multiple dimensionality reduction techniques, ranging from classical approaches to modern neural network-based methods. Based on empirical testing, these methods can be broadly categorised by their effectiveness in taxonomic separation:

### Highly Effective Methods

#### Principal Component Analysis (PCA) and Variants

- Standard PCA: Linear dimensionality reduction projecting data onto directions of maximum variance [1]
- PCA with SVD solver: Uses randomised Singular Value Decomposition for improved performance with large datasets ([Implementation](https://www.tutorialspoint.com/scikit_learn/scikit_learn_dimensionality_reduction_using_pca.htm)) [2]
- Kernel PCA: Non-linear variant using the kernel trick for capturing more complex relationships ([Implementation](https://www.tutorialspoint.com/scikit_learn/scikit_learn_dimensionality_reduction_using_pca.htm)) [3]

Key advantages:

- Computationally efficient
- Deterministic results (same output for same input)
- Well-understood mathematical properties
- Consistently reliable for taxonomic separation

#### Uniform Manifold Approximation and Projection (UMAP)

- Non-linear dimensionality reduction that preserves both local and global structure [4]
- Particularly effective at separating distinct taxonomic groups while maintaining relationships between similar sequences [5]

Key advantages:

- Better preservation of global structure compared to t-SNE
- Faster computation than t-SNE for large datasets
- Creates well-defined, compact clusters with natural shapes

#### Autoencoder with UMAP

- Neural network-based approach combining deep learning with UMAP visualisation [6]
- Particularly useful for complex datasets
- See [Autoencoder Documentation](kmers_autoencoder.md) for detailed implementation

### Moderately Effective Methods

#### t-Distributed Stochastic Neighbor Embedding (t-SNE)

- Non-linear technique emphasising local structure preservation [7]
- Fixed at 3 components due to library constraints
- May not preserve global structure as well as UMAP
- Significantly slower than PCA, especially for larger datasets
- Can take hours to compute for datasets with thousands of sequences

#### Isomap

- Non-linear dimensionality reduction using geodesic distances [8]
- Preserves global geometry
- Works well with data lying on a manifold
- Can be computationally intensive for large datasets
- Key parameters:
  - n_neighbors: Number of neighbors to consider (optimized automatically)
  - n_components: Number of dimensions in output (typically 2-3)
  - metric: Distance metric used (default: 'euclidean')

#### Random Trees Embedding

- Unsupervised transformation using randomised decision trees [9]
- Implementation based on [scikit-learn manifold learning examples](https://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html#sphx-glr-auto-examples-manifold-plot-lle-digits-py)
- Handles non-linear relationships
- Memory-efficient for large datasets

#### Non-Negative Matrix Factorization (NMF)

- Decomposes data into non-negative components [10]
- Implementation based on [scikit-learn NMF](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html)
- Natural for count data like k-mer frequencies
- May require parameter tuning

### Methods with Limited Testing

#### Locally Linear Embedding (LLE) Variants

- Standard LLE: Preserves local geometry by linear reconstruction [11]
- Modified and Hessian variants available
- Implementation based on [scikit-learn manifold learning examples](https://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html#sphx-glr-auto-examples-manifold-plot-lle-digits-py)
- Less reliable for taxonomic separation
- Sensitive to parameter choice
- Key parameters:
  - n_neighbors: Number of neighbors to consider (optimized automatically)
  - n_components: Number of dimensions in output
  - reg: Regularization constant
  - eigen_solver: Method for eigenvalue decomposition

#### Multidimensional Scaling (MDS)

- Preserves pairwise distances between points [12]
- Computationally intensive
- Limited performance in taxonomic separation

#### Spectral Embedding

- Based on eigendecomposition of similarity matrix [13]
- May not handle noise well
- Inconsistent results with genomic data
- Key parameters:
  - n_neighbors: Number of neighbors in adjacency matrix (optimized automatically)
  - n_components: Number of dimensions in output
  - affinity: Type of affinity matrix to construct

## Method Selection Guide

### Quick Start Recommendation

For most datasets, we recommend starting with **PCA**:

- Fast and reliable
- Consistently good separation
- Easy to interpret

### Method Selection Guide

1. **Simple Dataset or Initial Analysis**

   - Use PCA
   - Quick results
   - Well-understood properties

2. **Need for More Detailed Analysis**

   - Try UMAP
   - Better preservation of both local and global structure
   - More natural cluster shapes

3. **Complex Dataset or Unclear Separation**
   - Consider Autoencoder + UMAP
   - More computational power required
   - Can reveal subtle patterns

## Interpreting Results

### Component Selection

Different components may reveal different aspects of sequence relationships:

- Examine multiple component combinations (e.g., dim1 vs dim2, dim2 vs dim3)
- The most informative components may vary between datasets
- PC1 often correlates with GC content in k-mer based PCA

### Method Characteristics

**PCA**

- Reliable first choice for analysis
- Linear separation of taxonomic groups
- Quick to compute and interpret

**UMAP**

- Creates well-defined, compact clusters
- Preserves both local and global structure
- Clusters have natural, organic shapes

**Autoencoder + UMAP**

- Combines deep learning with effective visualisation
- Useful for complex datasets
- More computationally intensive

## Pipeline Configuration

### YAML Configuration Options

#### k-mer Length

- Specified using the `kmer_len` key in the ASCC input YAML file
- Integer value specifying k-mer size
- Default: 7
- Affects both computation time and memory usage
- Larger values provide more specific patterns but require more memory

#### Dimensionality Reduction Methods

The selection of kmer count dimensionality reduction methods that will be run is specified in with the 'dimensionality_reduction_methods' key in the ASCC input YAML file. The methods are provided in the YAML file as a list, e.g. `["pca", "umap", "t-sne"]`.

Available methods:

- `pca`
- `umap`
- `t-sne`
- `isomap`
- `lle_standard`
- `lle_hessian`
- `lle_modified`
- `mds`
- `se`
- `random_trees`
- `kernel_pca`
- `pca_svd`
- `nmf`
- `autoencoder_relu`
- `autoencoder_sigmoid`
- `autoencoder_linear`
- `autoencoder_tanh`
- `autoencoder_selu`

### Pipeline Integration

- Activated by including `kmers` in the `--include` parameter when running the ASCC pipeline
- Automatically handles file paths and integration with BlobToolKit
- Manages intermediate files and results

## `n_neighbors` Parameter Optimisation

### Purpose

Automatically finds optimal `n_neighbors` value for applicable methods:

- UMAP
- Isomap
- LLE (`standard`, `hessian`, `modified`)
- Spectral Embedding

### How it Works

- Tests multiple `n_neighbors` values (typically 5, 10, 15, 20, 30)
- For each value:
  - Performs dimensionality reduction
  - Applies HDBSCAN clustering
  - Calculates silhouette scores
- Selects `n_neighbors` value with best average clustering quality

### When to Use

- Enabled by default (can be disabled with `--skip_n_neighbors_optimisation` of the `kmer_count_dim_reduction.py` script)
- Particularly useful when:
  - Dealing with new types of sequence data
  - Dataset size differs significantly from previous analyses
  - Optimal neighbourhood size is unknown

### Output

- Generates optimisation results file (`method_name_n_neighbors_optimisation.txt`)
  - Contains tested `n_neighbors` values and their scores
  - Helps understand parameter selection process
- Logs best `n_neighbors` value in console output
- Uses the best `n_neighbors` value for generating the final dimensionality reduction embedding that is passed on to the downstream processes

## References

1. Jolliffe IT, Cadima J. Principal component analysis: a review and recent developments. Philos Trans A Math Phys Eng Sci. 2016
2. Halko N, Martinsson PG, Tropp JA. Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions. SIAM Review. 2011
3. Schölkopf B, Smola A, Müller KR. Nonlinear Component Analysis as a Kernel Eigenvalue Problem. Neural Computation. 1998
4. McInnes L, Healy J, Melville J. UMAP: Uniform Manifold Approximation and Projection for Dimension Reduction. ArXiv. 2018
5. McInnes L, Healy J, Saul N, Großberger L. UMAP: Uniform Manifold Approximation and Projection. Journal of Open Source Software. 2018
6. Hinton GE, Salakhutdinov RR. Reducing the Dimensionality of Data with Neural Networks. Science. 2006
7. van der Maaten L, Hinton G. Visualizing Data using t-SNE. Journal of Machine Learning Research. 2008
8. Tenenbaum JB, de Silva V, Langford JC. A Global Geometric Framework for Nonlinear Dimensionality Reduction. Science. 2000
9. Geurts P, Ernst D, Wehenkel L. Extremely Randomized Trees. Machine Learning. 2006
10. Lee DD, Seung HS. Learning the Parts of Objects by Non-negative Matrix Factorization. Nature. 1999
11. Roweis ST, Saul LK. Nonlinear Dimensionality Reduction by Locally Linear Embedding. Science. 2000
12. Borg I, Groenen PJF. Modern Multidimensional Scaling: Theory and Applications. Springer. 2005
13. Ng AY, Jordan MI, Weiss Y. On Spectral Clustering: Analysis and an Algorithm. NIPS. 2001
