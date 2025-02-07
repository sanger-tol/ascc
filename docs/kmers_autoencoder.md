# Variational Autoencoder Implementation for K-mer Analysis

## Overview

The pipeline implements a specialised Variational Autoencoder (VAE) for genomic sequence analysis, building upon the architecture described in Weber (2024) [1] for contamination detection in long-read sequencing data. While maintaining the core VAE principles from Weber's work, this implementation includes several modifications for improved scalability and performance, including dynamic layer sizing, memory-efficient processing, and enhanced visualisation capabilities through Uniform Manifold Approximation and Projection (UMAP) integration.

### Key Features

- Dynamic architecture that adapts to input data dimensions
- Emerging clusters in the latent space are likely to reflect taxonomic relationships, aiding downstream analysis
- Memory-efficient processing for large genomic datasets
- Integration with UMAP for improved visualisation
- Comprehensive quality monitoring and validation

## Architecture Details

### Input Layer and Preprocessing

- Accepts normalised k-mer frequency vectors (typically 16,384 dimensions for 7-mers)
- Implements Term Frequency-Inverse Document Frequency ([TF-IDF](https://en.wikipedia.org/wiki/Tf%E2%80%93idf)) transformation to emphasise distinctive k-mer patterns
- Includes [batch normalisation](https://en.wikipedia.org/wiki/Batch_normalization) [3] for stable training
- Uses length normalisation to account for varying sequence sizes

### Encoder Structure

- Dynamic sizing: intermediate layers automatically scale based on input dimensions
- Three main dense layers with decreasing dimensions (e.g., 2048 → 512 → 128)
- Residual connections to preserve important sequence features
- [Regularisation](<https://en.wikipedia.org/wiki/Regularization_(mathematics)>) through [dropout](<https://en.wikipedia.org/wiki/Dropout_(neural_networks)>) [4] (10%) and [L1 regularisation](<https://en.wikipedia.org/wiki/Lasso_(statistics)>) to prevent overfitting
- Batch normalisation after each major layer

### Latent Space

- Default 32 dimensions (adjustable via `n_components` parameter)
- Uses [reparameterisation trick](https://arxiv.org/abs/1312.6114) [5] with learned mean and variance
- Implements [Kullback-Leibler (KL) divergence](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence) warm-up to prevent posterior collapse
- The learned latent space is likely to reflect taxonomic relationships between sequences, aiding downstream analysis

### Decoder Structure

- Mirror of encoder architecture with transposed dimensions
- Additional regularisation to ensure non-negative outputs
- Residual connections matching encoder structure
- Final sigmoid activation to reconstruct normalized k-mer frequencies between 0-1

### Activation Modes

The autoencoder supports multiple activation functions for its dense layers:

- **ReLU** (Rectified Linear Unit): A widely used default activation, effective in preventing vanishing gradients
- **SeLU** (Scaled Exponential Linear Unit): Self-normalizing, suitable for deeper networks
- **Tanh**: Maps inputs to [-1,1], useful for symmetric distributions
- **Sigmoid**: Compresses outputs between [0,1], potentially limiting gradient updates
- **Linear**: No transformation, mainly used for debugging

Empirical testing has shown that ReLU is the recommended default choice, with SeLU also performing well. Running all activation modes is generally unnecessary; ReLU is sufficient for most use cases.

## Key Differences from Standard VAEs

### Data-specific Modifications

- Non-negativity constraints suited for k-mer counts
- Dynamic architecture sizing based on input k-mer dimensionality
- Standard mean squared error reconstruction loss with KL divergence weighted by adaptive β scheduling

### Training Optimisations

- Batch processing to handle large genome datasets
- Memory requirements scale with dataset size (default 4GB, user-adjustable)
- Modified KL divergence scheduling for stable training
- Custom [early stopping](https://en.wikipedia.org/wiki/Early_stopping) based on validation loss

### Genomic-specific Features

- Integration with taxonomic analysis workflows
- Preservation of important genomic signals
- Robust handling of varying sequence lengths

## Training Process

### Optimisation Strategy

- Uses an [adaptive learning rate](https://en.wikipedia.org/wiki/Learning_rate#Adaptive_learning_rate) strategy with [Adam optimizer](https://arxiv.org/abs/1412.6980) [6]
- Implements early stopping to prevent overfitting
- Monitors reconstruction quality and latent space structure
- Automatically adjusts KL divergence weight

### Parameters and Tuning

- Latent space dimensions adjustable based on dataset complexity
- Batch size and epoch count adapt to dataset size
- Automatic prevention of common VAE issues like posterior collapse
- Comprehensive monitoring tools for performance tracking

## Output Files and Visualisation

### Training History Plot (`vae_training_history.png`)

Four key training metrics:

- Total Loss: Combined reconstruction and KL divergence loss
- Reconstruction Loss: Measures accuracy of sequence reconstruction
- KL Divergence Loss: Measures latent space distribution
- Learning Rate: Shows adaptive learning rate behaviour

### Embedding Visualisations

Generated for both raw VAE embeddings and UMAP transformations:

- Main scatter plot of the 2D embedding
- Density view showing point concentrations
- Distribution plots for first two components

### Training Log (`training_log.csv`)

Detailed per-epoch metrics including:

- Total loss
- Reconstruction loss
- KL divergence loss
- Validation metrics
- Learning rate

### Quality Metrics

Files containing embedding quality measures:

- Trustworthiness: Measures how well local neighborhoods are preserved in the embedding
  - Values range from 0 to 1
  - Values >0.9: Excellent preservation of local structure
  - Values 0.7-0.9: Good preservation, suitable for most analyses
  - Values <0.7: Poor preservation, consider adjusting parameters or using a different method
- Continuity: Measures whether points that are close in the original space remain close in the embedding
  - Values range from 0 to 1
  - Values >0.9: Strong preservation of original relationships
  - Values 0.7-0.9: Acceptable preservation of relationships
  - Values <0.7: Significant distortion, may not be reliable for analysis
- Both metrics should ideally be >0.7 for reliable embeddings
- If either metric is consistently low, consider:
  - Adjusting the number of components
  - Using a different activation function
  - Trying a different dimensionality reduction method

### Coordinate Files

- `autoencoder_[activation mode]_raw_coordinates.csv`: Raw embeddings from the autoencoder, where the file name reflects the activation mode used
- `kmers_dim_reduction_embeddings.csv`: Final embeddings combining VAE latent space with UMAP projection, optimized for species separation. UMAP optimization occurs when `--skip_n_neighbors_optimisation` isn't used (default).

## Performance Considerations

### Memory Management

- Processes large genomic datasets through batch processing
- Memory usage increases with dataset size and k-mer length
- Default memory limit of 4GB (adjustable via configuration)
- Automatically scales intermediate layer dimensions
- Uses data normalisation strategies optimised for k-mer frequency data
- Maintains numerical stability through appropriate regularisation
- For datasets with <100 sequences, autoencoder components automatically reduce to sample_count-1
- Warning: Large datasets may require significant memory resources

### Integration with UMAP

- Combines VAE's dimensionality reduction with UMAP's visualisation
- Preserves both local and global structure
- Enables interactive exploration of sequence relationships
- Provides options for different visualisation schemes

## Interpreting Results

While traditional GC-coverage plots provide a basic view of sequence composition, k-mer-based embeddings often reveal more nuanced taxonomic relationships. The first principal component (PC1) frequently correlates with GC content, but additional components and non-linear transformations can expose more subtle patterns in sequence composition that aid in taxonomic classification.

## Applications

### Primary Use Cases

- Separating sequences from different organisms in mixed samples
- Detecting contamination in sequencing data
- Exploring relationships between sequences
- Identifying unusual or outlier sequences

### Limitations and Considerations

- Very similar sequences might not be completely separated
- Processing time increases with dataset size and k-mer length
- Memory requirements scale with input dimensions
- Training time can be significant for large datasets

## References

1. Weber, Claudia C. "Disentangling cobionts and contamination in long-read genomic data using sequence composition." G3 Genes|Genomes|Genetics, 2024.
2. Hinton GE, Salakhutdinov RR. "Reducing the Dimensionality of Data with Neural Networks." Science, 2006.
3. Ioffe S, Szegedy C. "Batch Normalization: Accelerating Deep Network Training by Reducing Internal Covariate Shift." arXiv:1502.03167, 2015.
4. Srivastava N, Hinton G, Krizhevsky A, Sutskever I, Salakhutdinov R. "Dropout: A Simple Way to Prevent Neural Networks from Overfitting." Journal of Machine Learning Research, 2014.
5. Kingma DP, Welling M. "Auto-Encoding Variational Bayes." arXiv:1312.6114, 2013.
6. Kingma DP, Ba J. "Adam: A Method for Stochastic Optimization." arXiv:1412.6980, 2014.
