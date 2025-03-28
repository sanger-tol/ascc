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

### Data Normalisation

Before applying dimensionality reduction methods, the k-mer count data is normalised to a range of 0 to 1. This normalisation step:

- Ensures all features contribute equally to the analysis
- Improves numerical stability during computation
- Makes the results more comparable across different sequences
- Helps prevent any single feature from dominating the analysis

## Validation Dataset

To demonstrate and evaluate the various dimensionality reduction methods, we use a synthetic dataset derived from complete genome assemblies. This dataset was specifically constructed to represent a realistic scenario of mixed genomic content while maintaining reliable taxonomic labels for validation purposes.

### Dataset Construction

The dataset was created by collecting complete genome sequences from diverse organisms:

- **Chlamydomonas reinhardtii**:
  - Nuclear genome (CM023806.1)
  - Chloroplast (NC_005353.1)
  - Mitochondrion (NC_001638.1)
- **Aspergillus fumigatus**:
  - Nuclear genome (NC_007194.1)
  - Mitochondrion (NC_017016.1)
- **Pseudooceanicola atlanticus** (CP051248.1)
- **Streptococcus agalactiae** (CP051848.1)
- **Brevibacterium siliguriense** (NZ_LT629766.1)
- **Halothece sp. PCC 7418** (NC_019779.1)
- **Methanofollis liminatans** (NZ_CM001555.1)

The sequences were fragmented with:

- Maximum length of 50 kb
- Resulting in 3,177 sequence fragments
- Total assembly size of 158 Mb
- Average fragment length of ~49.7 kb

### Dataset Characteristics

The dataset contains a taxonomically diverse mix of sequences:

- Eukaryotic (algal and fungal)
- Bacterial
- Archaeal
- Organellar (chloroplast and mitochondrial)

Key features:

- Represents multiple distinct taxonomic groups in a single dataset
- Includes both nuclear and organellar DNA from the same species
- Comprises sequences of uniform length (mostly 50 kb, with shorter fragments at contig ends)

### Rationale for Synthetic Dataset

This synthetic dataset was chosen over real environmental samples for several key reasons:

1. **Reliable Sequence Classifications**:

   - All sequences are publicly available
   - Pre-existing, high-confidence taxonomic classifications from NCBI

2. **Validation Independence**:

   - External classifications prevent circular validation when testing taxonomic separation methods

3. **Fragment Size Control**:

   - Uniform fragmentation helps evaluate method performance across consistent sequence lengths
   - Short sequences (≤50 kb) represent the most difficult cases for taxonomic classification

4. **Statistical Power**:
   - The fragmentation process provides sufficient sequence count (>3,000)
   - Enables robust evaluation of dimensionality reduction methods

The dataset's characteristics make it particularly suitable for evaluating taxonomic separation methods, as it combines the complexity of multi-species genomic content with reliable reference classifications for validation.

## Available Methods

The tool implements multiple dimensionality reduction techniques, ranging from classical approaches to modern neural network-based methods. Based on empirical testing with our validation dataset, these methods can be broadly categorised by their effectiveness in taxonomic separation. Example visualisations for each method are available in the `kmers_dim_reduction_btk_figures` directory.

### Highly Effective Methods

#### Principal Component Analysis (PCA) and Variants

- Standard PCA: Linear dimensionality reduction projecting data onto directions of maximum variance [1] ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.pca.png))
- PCA with SVD solver: Uses randomised Singular Value Decomposition for improved performance with large datasets ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.pca_with_svd_solver.png)) [2]
- Kernel PCA: Non-linear variant using the kernel trick for capturing more complex relationships ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.kernel_pca.png)) [3]

Key advantages:

- Computationally efficient
- Deterministic results (same output for same input)
- Well-understood mathematical properties
- Consistently reliable for taxonomic separation

#### Uniform Manifold Approximation and Projection (UMAP)

- Non-linear dimensionality reduction that preserves both local and global structure [4] ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.umap.png))
- Particularly effective at separating distinct taxonomic groups while maintaining relationships between similar sequences [5]

Key advantages:

- Better preservation of global structure compared to t-SNE
- Faster computation than t-SNE for large datasets
- Creates well-defined, compact clusters with natural shapes

#### Autoencoder with UMAP

- Neural network-based approach combining deep learning with UMAP visualisation [6]
- Particularly useful for complex datasets
- Multiple activation functions available:
  - ReLU ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.autoencoder_relu_umap.png))
  - SELU ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.autoencoder_selu_umap.png))
  - Tanh ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.autoencoder_tanh_umap.png))
  - Linear (BlobToolKit visualisation unavailable due to UI limitation with long text overlapping control buttons)
- See [Autoencoder Documentation](kmers_autoencoder.md) for detailed implementation

### Moderately Effective Methods

#### t-Distributed Stochastic Neighbor Embedding (t-SNE)

- Non-linear technique emphasising local structure preservation [7] ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.t-SNE.png))
- Fixed at 3 components due to library constraints
- May not preserve global structure as well as UMAP
- Significantly slower than PCA, especially for larger datasets
- Can take hours to compute for datasets with thousands of sequences

#### Isomap

- Non-linear dimensionality reduction using geodesic distances [8] ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.isomap.png))
- Preserves global geometry
- Works well with data lying on a manifold
- Can be computationally intensive for large datasets
- Key parameters:
  - `n_neighbors`: Number of neighbours to consider (optimised automatically)
  - n_components: Number of dimensions in output (typically 2-3)
  - metric: Distance metric used (default: 'euclidean')

#### Random Trees Embedding

- Unsupervised transformation using randomised decision trees [9] ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.random_trees.png))
- Implementation based on [scikit-learn manifold learning examples](https://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html#sphx-glr-auto-examples-manifold-plot-lle-digits-py)
- Handles non-linear relationships
- Memory-efficient for large datasets

#### Non-Negative Matrix Factorization (NMF)

- Decomposes data into non-negative components [10] ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.non_negative_matrix_factorisation.png))
- Implementation based on [scikit-learn NMF](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html)
- Natural for count data like k-mer frequencies
- May require parameter tuning

### Methods with Limited Testing

#### Locally Linear Embedding (LLE) Variants

- Standard LLE: Preserves local geometry by linear reconstruction [11] ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.LLE_standard.png))
- Hessian variant available ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.LLE_hessian.png))
- Implementation based on [scikit-learn manifold learning examples](https://scikit-learn.org/stable/auto_examples/manifold/plot_lle_digits.html#sphx-glr-auto-examples-manifold-plot-lle-digits-py)
- Less reliable for taxonomic separation
- Sensitive to parameter choice
- Key parameters:
  - `n_neighbors`: Number of neighbours to consider (optimised automatically)
  - n_components: Number of dimensions in output
  - reg: Regularisation constant
  - eigen_solver: Method for eigenvalue decomposition

#### Multidimensional Scaling (MDS)

- Preserves pairwise distances between points [12] ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.MDS.png))
- Computationally intensive
- Limited performance in taxonomic separation

#### Spectral Embedding

- Based on eigendecomposition of similarity matrix [13] ([example plot](kmers_dim_reduction_btk_figures/btk_datasets_CBD.blob.circle.spectral_embedding.png))
- May not handle noise well
- Inconsistent results with genomic data
- Key parameters:
  - `n_neighbors`: Number of neighbours in adjacency matrix (optimised automatically)
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

## BlobToolKit Visualisation Results

This section presents BlobToolKit visualisations of the Chlamydomonas reinhardtii test dataset using various dimensionality reduction methods. The methods are organised by their effectiveness in separating taxonomic groups in this specific dataset.

### High-performing Methods

These methods demonstrated good separation of taxonomic groups in the test dataset, creating well-defined clusters that align with known taxonomic classifications.

#### PCA and Variants

![PCA](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.pca.png)

Standard PCA provides clear separation of taxonomic groups in the Chlamydomonas dataset, with distinct clusters for the different organisms and organelles.

![PCA with SVD Solver](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.pca_with_svd_solver.png)

PCA with SVD solver produces results similar to standard PCA but can be more efficient for large-scale genomic data.

![Kernel PCA](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.kernel_pca.png)

Kernel PCA captures more complex relationships between sequences, resulting in well-defined clusters for the different taxonomic groups.

#### UMAP

![UMAP](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.umap.png)

UMAP creates compact, well-separated clusters for the different taxonomic groups in the dataset, demonstrating its effectiveness for taxonomic separation.

#### Autoencoder with UMAP

![Autoencoder Linear + UMAP](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle_autoencoder_linear_umap.png)

The combination of a linear autoencoder with UMAP produces good separation of taxonomic groups.

![Autoencoder ReLU + UMAP](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.autoencoder_relu_umap.png)

Autoencoder with ReLU activation followed by UMAP creates distinct clusters for different taxonomic groups.

![Autoencoder SELU + UMAP](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.autoencoder_selu_umap.png)

Autoencoder with SELU activation followed by UMAP also performs well in separating taxonomic groups.

![Autoencoder Tanh + UMAP](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.autoencoder_tanh_umap.png)

Autoencoder with Tanh activation followed by UMAP demonstrates good separation of taxonomic groups.

#### Isomap

![Isomap](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.isomap.png)

Isomap effectively separates the taxonomic groups in the dataset, preserving global geometry.

#### Spectral Embedding

![Spectral Embedding](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.spectral_embedding.png)

Spectral Embedding performs well in the Chlamydomonas dataset, showing clear separation between taxonomic groups despite being categorised as a method with limited testing in general evaluations.

### Intermediate-performing Methods

These methods show some separation between taxonomic groups, but the boundaries are less clear and clusters may overlap.

#### t-SNE

![t-SNE](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.t-SNE.png)

t-SNE shows separation between taxonomic groups, but the boundaries are somewhat fuzzy and clusters overlap.

#### Locally Linear Embedding (Standard)

![LLE Standard](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.LLE_standard.png)

Standard LLE provides separates some of the species well but others clump together.

#### Random Trees Embedding

![Random Trees](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.random_trees.png)

Random Trees Embedding has separated the eukaryotic species (Chlamydomonas and Aspergillus) well. The bacteria are separated from the eukaryotic species but the clusters of individual bacterial species have some overlap with one another.

### Low-performing Methods

These methods did not effectively separate taxonomic groups in the test dataset.

#### Autoencoder with Sigmoid

![Autoencoder Sigmoid](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.autoencoder_sigmoid.png)

Autoencoder with Sigmoid activation does not effectively separate taxonomic groups in the dataset.

#### Non-Negative Matrix Factorization (NNMF)

![NNMF](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.nnmf.png)

NNMF fails to create distinct clusters for different taxonomic groups in this dataset.

#### Multidimensional Scaling (MDS)

![MDS](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.MDS.png)

MDS does not effectively separate taxonomic groups in the dataset.

#### Locally Linear Embedding (Hessian)

![LLE Hessian](kmers_dim_reduction_btk_plots/btk_datasets_CBD.blob.circle.LLE_hessian.png)

Hessian LLE performs poorly in separating taxonomic groups in this dataset.

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
