# General guidelines

## Course goal

This course walks you through a complete metagenomics analysis across **five days**. On **Day 1** you will learn the core concepts of metagenomics, study-design considerations and sequencing technologies, and how to assess the quality of raw sequencing reads (both short, Illumina-style reads and long, PacBio-style reads) with FastQC. On **Day 2** we cover targeted marker-gene sequencing: after the theory (marker gene selection, primer design, PCR biases, OTUs vs ASVs), you will process 16S amplicon (eDNA) reads with DADA2 (filtering, error correction, merging, chimera removal, taxonomy assignment), followed by a downstream analysis with PhyloSeq (rarefaction, α- and β-diversity, barplots, geographic localisation and a LEfSe analysis). On **Day 3** we move to shotgun (whole-genome) metagenomics: how depth and coverage shape the analysis, k-mer vs alignment-based read classification, and a Nextflow-based pipeline that performs host removal with Bowtie2, taxonomic classification with Kraken2/Bracken and visualization with Krona plots. You will then compare samples with Pavian and build a Phyloseq object from the pipeline output to explore α- and β-diversity, ordination and network analyses. On **Day 4** we deal with assembly and binning: principles of _de novo_ metagenome assembly, the challenges of community complexity and strain microdiversity, and the binning strategies (MetaBAT2, MaxBin2 and friends) used by nf-core/mag to recover Metagenome-Assembled Genomes. On **Day 5** we close the analysis arc with the integral analysis of the recovered MAGs — quality assessment with BUSCO/CheckM2 and the BIgMAG dashboard, metabolic potential with KEGG Decoder, functional annotation with COGclassifier, dbCAN3 and antiSMASH, genome visualization with Proksee and a pangenome comparison with Roary/Phandango — before zooming out to FAIR data and reproducible analyses with workflow managers and containers.

By the end of the course, you will have run a complete metagenomics workflow from raw reads to biological interpretation, and you will be able to interpret the outputs at each stage: FastQC reports, DADA2 denoising outputs, Kraken2/Bracken reports, Krona plots, Pavian comparisons, Phyloseq figures and MAG quality/functional annotation outputs.

## Software

Almost all the software needed in this course is either:

* Already installed in the GitHub Codespaces environment.
* Installed on the fly with Conda during the exercises.
* Available as R/Bioconductor packages that you install once beforehand.

The `R` scripts in the practicals typically require the Phyloseq ecosystem. You can install all the packages used throughout the course with the script below:

[:fontawesome-solid-file-code: Download the required packages script](../assets/code/required_packages.R){: .md-button }

??? full-code "required_packages.R"
    ```r title="required_packages.R" linenums="1"
    # these commands should be run once
    if (!require("BiocManager", quietly = TRUE))

      install.packages("BiocManager")
    BiocManager::install(version = "3.23")

    install.packages("devtools") install.packages("sf")
    install.packages("magick")
    install.packages("mapview")

    BiocManager::install("dada2")
    BiocManager::install("phyloseq")
    BiocManager::install("lefser")
    BiocManager::install("Biostrings")
    BiocManager::install("ggplot2")
    BiocManager::install("tidyr")
    BiocManager::install("dplyr")
    BiocManager::install("magrittr")
    BiocManager::install("vegan")
    BiocManager::install("gridExtra")
    #BiocManager::install("microbiomeMarker")
    devtools::install_github("zachcp/phylogeo")

    # These might require some recompilation and ask for additional packages,
    # please answer Yes to the questions.

    # You can then check for each library if they load with these commands:

    library(dada2, quietly = TRUE)
    library(phyloseq, quietly = TRUE)
    library(Biostrings, quietly = TRUE)
    library(ggplot2, quietly = TRUE)
    library(magrittr, quietly = TRUE)
    library(vegan, quietly = TRUE)
    library(gridExtra, quietly = TRUE)
    library(phylogeo, quietly = TRUE)

    # problem if microbiomeMarker does not exist for Bioc 3.23...

    install.packages("remotes")
    remotes::install_github("yiluheihei/microbiomeMarker")

    library(microbiomeMarker, quietly = TRUE)
    ```

!!! tip "Hands-on environment"
    All the exercises can be run either locally (RStudio/VS Code + R) or on a GitHub Codespaces instance. We provide the data files (`test_sheet.csv`, `merged.biom`) and the Krona reports as static assets on this website so that you can reproduce the analysis without re-running the whole pipeline.

## Website colour code explanation

We tried to use a colour code throughout the website to make the different pieces of information easily distinguishable. Here's a quick summary about the colour blocks you will encounter:

!!! info "This is a supplementary piece of information"

!!! tip "This is a tip to help you solve an exercise"

!!! success "This is the answer to an exercise"

!!! warning "This is a warning about a potential problem"

!!! bug "This is an explanation about a common bug/error"

!!! full-code "This is a full code file"

!!! abstract "These are directory contents"

## Exercises

Each question provides a background explanation, a description of the task at hand and additional details when required.

!!! tip "Hints for challenging questions"
    For the most challenging questions, hints will be provided. However, you should first try to solve the problem without them!

## Answers

Do not hesitate to modify and overwrite your code from previous answers as difficulty is incremental. The questions are designed to incite you to build your answers upon the previous ones.

If something is not clear at any point, please raise your hand and we will do our best to answer your questions!

## Learning outcomes

**After having completed this course you will be able to:**

* Explain the core concepts of metagenomics, design a study (replication, controls, sampling, sequencing strategy) and distinguish targeted from shotgun approaches.
* Run FastQC on short and long reads and interpret its main report modules.
* Process targeted 16S amplicon reads with DADA2 and analyse the recovered ASVs with PhyloSeq (rarefaction, α- and β-diversity, geographic mapping and LEfSe).
* Explain how sequencing depth and coverage shape a shotgun metagenomics analysis, and how k-mer- and alignment-based classifiers work.
* Run a reference-based metagenomics characterization pipeline (Bowtie2 → Kraken2/Bracken → Krona).
* Compare taxa across samples interactively with Krona plots and Pavian.
* Build a Phyloseq object from a BIOM file and sample metadata.
* Estimate α- and β-diversity and run ordination and network analyses.
* Explain the principles of _de novo_ metagenome assembly and evaluate assembly quality.
* Describe binning strategies and their role in recovering MAGs.
* Assess the quality and taxonomy of MAGs with BUSCO/CheckM2 and the BIgMAG dashboard.
* Infer the metabolic potential of MAGs with KEGG Decoder.
* Annotate functional categories (COG, CAZymes, biosynthetic gene clusters), visualize genomes with Proksee and run a pangenome analysis.
* Apply the FAIR principles and use workflow managers/containers to make metagenomics analyses reproducible.