# General guidelines

## Course goal

<p style='text-align: justify;'>This course walks you through a complete metagenomics analysis. During the **first session** you will learn how to assess the quality of raw sequencing reads (both short, Illumina-style reads and long, PacBio-style reads) with FastQC and how to interpret the different modules of the FastQC report. During the **second session**, you will process targeted 16S amplicon (eDNA) reads with DADA2: filtering, error correction, merging, chimera removal and taxonomy assignment, followed by a downstream analysis with PhyloSeq (rarefaction, α- and β-diversity, barplots, geographic localisation and a LEfSe analysis). During the **third session**, you will run a Nextflow-based pipeline that performs host removal with Bowtie2, taxonomic classification with Kraken2/Bracken, and visualization with Krona plots. You will then compare samples with Pavian and move to a downstream compositional analysis: building a Phyloseq object from the pipeline output and exploring α- and β-diversity, ordination and network analyses. During the **fourth session**, you will perform an integral analysis of the Metagenome-Assembled Genomes (MAGs) recovered with nf-core/mag: overall quality assessment with the BIgMAG dashboard, inference of the metabolic potential with KEGG Decoder, functional annotation with COGclassifier, dbCAN3 and antiSMASH, genome visualization with Proksee and a pangenome comparison with Roary/Phandango.</p>

<p style='text-align: justify;'>By the end of the course, you will have run a complete metagenomics workflow from raw reads to biological interpretation, and you will be able to interpret the outputs at each stage: FastQC reports, DADA2 denoising outputs, Kraken2/Bracken reports, Krona plots, Pavian comparisons and Phyloseq figures.</p>

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

* Run FastQC on short and long reads and interpret its main report modules.
* Process targeted 16S amplicon reads with DADA2 and analyse the recovered ASVs with PhyloSeq (rarefaction, α- and β-diversity, geographic mapping and LEfSe).
* Run a reference-based metagenomics characterization pipeline (Bowtie2 → Kraken2/Bracken → Krona).
* Compare taxa across samples interactively with Krona plots and Pavian.
* Build a Phyloseq object from a BIOM file and sample metadata.
* Estimate α- and β-diversity and run ordination and network analyses.
* Assess the quality and taxonomy of MAGs with the BIgMAG dashboard.
* Infer the metabolic potential of MAGs with KEGG Decoder.
* Annotate functional categories (COG, CAZymes, biosynthetic gene clusters), visualize genomes with Proksee and run a pangenome analysis.