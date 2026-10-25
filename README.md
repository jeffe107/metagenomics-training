# Metagenomic Data Analysis of Microbial Communities (MGMC01)

![Static Badge](https://img.shields.io/badge/release%20date-october%202026-%230DC09D)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

- Website: [https://jeffe107.github.io/metagenomics-training/](https://jeffe107.github.io/metagenomics-training/)
- Issues: [https://github.com/jeffe107/metagenomics-training/issues](https://github.com/jeffe107/metagenomics-training/issues)

## Authors

- [Jeferyd Yepes-García](https://jeferydyepes.com/)

## Course content

- **Day 1 – Introduction and quality control**: Metagenomics concepts, study design, sequencing technologies (short and long reads). Quality control with FastQC for Illumina, PacBio and metagenomic datasets.
- **Day 2 – Targeted sequencing analysis**: Amplicon sequencing (marker genes, primer design, PCR biases, OTUs vs ASVs). Processing of 16S amplicon (eDNA) reads with DADA2 (filtering, error correction, merging, chimera removal, taxonomy assignment). Downstream analysis with Phyloseq (rarefaction, α- and β-diversity, barplots, geographic mapping, LEfSe).
- **Day 3 – Shotgun metagenomics and downstream analysis**: Principles of whole-genome shotgun sequencing (depth, coverage, classification vs profiling). Reference-based classification with a Nextflow pipeline (Bowtie2 host removal, Kraken2, Bracken), interactive Krona plots, Pavian exploration, and downstream Phyloseq analysis (α- and β-diversity, ordination, networks).
- **Day 4 – Metagenome assembly and binning**: _De novo_ metagenome assembly, quality metrics, binning (MetaBAT2, MaxBin2), and exploration of MAGs with nf-core/mag.
- **Day 5 – MAG annotation and FAIR data**: MAG quality assessment (BUSCO/CheckM2, BIgMAG), metabolic potential (KEGG Decoder), functional annotation (COG, CAZymes, biosynthetic gene clusters), genome visualization (Proksee), pangenome analysis, and FAIR data/reproducibility with workflow managers and containers.