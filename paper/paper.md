---
title: 'ASCC: a Nextflow pipeline for identifying cobionts and contaminants in genome assemblies'
tags:
  - nextflow
  - genomics
  - decontamination
  - pipeline
  - biodiversity
authors:
  - name: Eerik Aunin
    orcid: 0000-0001-8385-2636
    equal-contrib: true
    affiliation: 1
  - name: Damon-Lee B Pointon
    orcid: 0000-0003-2949-6719
    equal-contrib: true
    affiliation: 1
  - name: James Torrance
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: 1
  - name: Ying Sims
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Will Eagles
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Matthieu Muffato
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Shane A. McCarthy
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: The Sanger Institute, Wellcome Genome Campus, UK
   index: 1
date: 01 January 2025
bibliography: paper.bib

---

# Summary

When producing genome assemblies, a routine procedure is to check the assembled sequences for cobionts and contaminants. Numerous existing tools can be used for this purpose, however, each tool can produce spurious results. When dealing with assemblies from samples with a complex composition of species, it can be useful to run multiple tools for sequence identification in parallel and then compare their results. We have put together a Nextflow pipeline called Assembly Screen for Cobionts and Contaminants (ASCC) that collects a selection of taxonomic identification and decontamination tools in one pipeline. The pipeline contains Tiara, read mapping, BLAST and Kraken2 against the NCBI nt database, Diamond BLASTX against the NCBI nr and Uniprot databases, kmer counting and dimensionality reduction, BlobToolKit, FCS-GX, FCS-adaptor, VecScreen, a PacBio barcodes check and BLAST for detecting organellar sequences. Running each component of the pipeline is optional. The results of a run are collected as a BlobToolKit dataset and CSV tables. This pipeline is regularly used in the Tree of Life programme of the Wellcome Sanger Institute, running a large-scale assembly production pipeline.

# Statement of need

When generating genome assemblies, a major issue is the taxonomic identification of the assembled sequences. Assemblies of environmental samples can contain sequences from hundreds of species. Contaminant sequences can also appear in assemblies of samples expected to contain only one target species. Various existing tools and databases can be used for the identification of cobionts and contaminants in genome assemblies. All the individual tools have their strengths and weaknesses and can produce sequence classifications that disagree with results from other tools. Running multiple contaminant detection tools in parallel on the same assembly and comparing the results can help to determine which classifications are likely spurious and which are likely accurate. Moreover, the results from each tool have a different profile of missed identifications, so combining the output of multiple tools can reduce the number of unidentified sequences.

The ASCC (Assembly Screen for Cobionts and Contaminants) pipeline is designed to runrun multiple decontamination tools and then merge their results. One of the main outputs of the pipeline is a CSV (comma-separated values) file where the rows correspond to sequences in the assembly and the columns contain the classification results from the constituent tools.
Most of the existing software tools for identifying sequences in assemblies or reads are focused on prokaryotic species only. These include CheckM (Parks et al. 2015), Anvi’o (Eren et al. 2021), CLARK (Ounit et al. 2015), and GUNC (Orakov et al. 2021). Bacterial metagenomics tools are often designed for short-read data. ASCC differs from bacterial metagenomics tools, as its main use cases involve long-read assemblies of eukaryotic target species that may or may not contain bacterial sequences.
All individual components of the ASCC pipeline are optional when setting up a run. This allows for highly customised runs according to the assembly complexity and end-user needs, e.g. enabling lightweight runs where speed is essential. ASCC was developed using Nextflow (Di Tommaso et al. 2017) and follows the nf-core (Ewels et al. 2020) standards. This ensures adherence to pipeline development best practices and enables running the pipeline on a variety of compute platforms. While nf-core has existing pipelines for identifying species in reads (nf-core/mag (Krakau et al. 2022), nf-core/taxprofiler (Yates et al. 2023) and nf-core/detaxizer (Seidel et al. 2024)), there is currently no official nf-core pipeline for the decontamination of assemblies. The Nextflow decontamination pipeline CLEAN (Lataretu et al. 2023) only uses read mapping for identification of contaminants.
Some existing tools are meant for decontaminating genome assemblies but appear to be no longer maintained. An established and maintained tool used for identifying sequences in genome assemblies is BlobToolKit (Challis et al. 2020). One part of BlobToolKit is the BlobToolKit Pipeline (Sims et al. 2025), which uses BUSCO (Simão et al. 2015), Diamond (Buchfink, Reuter, and Drost 2021) and BLAST (Sayers et al. 2022) to identify the taxonomic origin of sequences. Another part of the BlobToolKit ecosystem is the web app BlobToolKit Viewer, which uses the JSON output from the BlobToolKit Pipeline's to produce an interactive visualisation. ASCC expands upon the BlobToolKit Pipeline by offering additional modules and therefore additional annotations to the output.

ASCC is currently used in the Tree of Life programme at the Wellcome Sanger Institute for decontamination of all produced assemblies, ranging from those assuming a single target species, to those generated from cobiont systems.


# Pipeline components

An overview of the components of the pipeline is shown in Figure 1.

## Inputs
The only essential input is an assembly FASTA file, multiple assemblies can be provided where there are primary, haplotype, mitochondrial and plastid assemblies (the requirement is that the assemblies belong to the same sample). The pipeline can take reads (both long and short reads) as an input and map them with minimap2 (Li 2018) to determine coverage, these reads are also used in the BlobToolKit pipeline, for more indepth analysis. ASCC also accepts organellar FASTA files as input to use in the detection of organellar contamination in a nuclear genome assembly (using BLAST) as well as some basic analysis to identify any contamination in the organellar genome. Version 1 of the ASCC pipeline will take a `params.yaml` detailing the assembly information, database locations and various software paramters.

## Data processing
Below are brief descriptions of the components of ASCC, as shown in Figure 1. FCS-GX (Astashyn et al. 2024) is NCBI's software for detecting contaminant species in assemblies using a cross-species genome aligner. FCS-adaptor (Astashyn et al. 2023) and VecScreen (NCBI 2001) are tools from NCBI for detecting adapter contamination. VecScreen is older than FCS-adaptor but its advantage over FCS-adaptor is that it allows the user to specify a custom database of adapter sequences. ASCC can use BLAST to detect contamination of the assembly with user-provided barcode sequences (a collection of common PacBio barcode sequences are supplied with the pipeline). Tiara (Karlicki, Antonowicz, and Karnkowska 2022) uses a neural network for domain-level taxonomic classification of sequences and the detection of organellar sequences.
ASCC can run the BlobToolKit pipeline. It can also run BLAST against the NCBI nt database and Diamond BLASTX against the NCBI nr and Uniprot (UniProt Consortium 2023) databases separately from the BlobToolKit Pipeline. With default settings, runs using ASCC's own BLAST and Diamond processes are faster and more lightweight compared to the BlobToolKit Pipeline. Besides BLAST against the NCBI nt database, ASCC can run Kraken (Wood, Lu, and Langmead 2019) with the nt database. Kraken uses kmer matches to identify sequences and is faster than BLAST. ASCC also has a component for kmer counting and dimensionality reduction of the kmer counts (using kcounter, Camargo. 2020), scikit-learn (Pedregosa et al. 2011) and Tensorflow (Abadi et al. 2015).

## Outputs
The outputs depend on which of sub-workflows of the pipeline were run, by default, this is set as ALL subworkflows. The key results from the pipeline are gathered as a CSV table and can be converted into a BlobToolKit dataset. BlobToolKit allows the user to incorporate custom variables into BlobToolKit datasets and then display them in BlobToolKit Viewer. ASCC makes use of this, creating BlobToolKit datasets that contain results from tools that are not a part of the BlobToolKit Pipeline. Running the BlobToolKit Pipeline is not required for making a BlobToolKit dataset with ASCC: it is possible to make a BlobToolKit dataset just using results from other components of ASCC.
If the required inputs are present, the pipeline can try to derive a phylum-level consensus classification from the outputs of multiple tools and report the estimated coverage per phylum. The adapter and organellar contamination reports are collected as plain text files.

# Availability and usage
ASCC and instructions for running it are available on GitHub at https://github.com/sanger-tol/ascc.

# Citations

- `@FCSgx` -> "(Astashyn et al. 2024)"
- `@FCSadaptor` -> "(Astashyn et al. 2023)"
- `@Buchfink_Reuter_Drost_2021` -> "(Buchfink, Reuter, and Drost. 2021)"
- `@Camargo` -> "Camargo. 2020"
- `@blobtoolkit` -> "(Challis et al. 2020)"
- `@btk_nxf` -> "(Sims et al. 2025)"
- `@Nextflow` -> "(Di Tommaso et al. 2017)"
- `@anvio` -> "(Eren et al. 2021)"
- `@nf-core` -> "(Ewels et al. 2020)"
- `@tiara` -> "(Karlicki, Antonowicz, and Karnkowska. 2022)"
- `@nf-core-mag` -> "(Krakau et al. 2022)"
- `@CLEAN` -> "(Lataretu et al. 2023)"
- `@minimap` -> "(Li. 2018)"
- `@tensorflow` -> "(Abadi et al. 2015)"
- `@vecscreen` -> "(NCBI. 2001)"
- `@gunc` -> "(Orakov et al. 2021)"
- `@clark` -> "(Ounit et al. 2015)"
- `@checkm` -> "(Parks et al. 2015)"
- `@scikit` -> "(Pedregosa et al. 2011)"
- `@ncbi_db` -> "(Sayers et al. 2022)"
- `@detaxiser` -> "(Seidel et al. 2024)"
- `@busco` -> "(Simão et al. 2015)"
- `@taxprofiler` -> "(Yates et al. 2023)"
- `@kraken2` -> "(Wood, Lu, and Langmead 2019)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

The creation of the ASCC pipeline was made possible by the Aquatic Symbiosis Genomics project, funded by the Moore Foundation. We thank Noah Gettle, Ksenia Krasheninnikova, Michael Paulini and Camilla Santos for testing the pipeline and reporting back issues encountered with it. We thank Kerstin Howe for her comments on the manuscript.

# References
