# Day 5 - FAIR principles, reproducibility and Project

## Learning outcomes

**After having completed this chapter you will be able to:**

* Identify different methodologies to enable reproducible research.
* Understand the requirements that every project must fulfill in order to be FAIR.
* Recognize the benefits of making projects Open Access.
* Describe the main steps to follow in order to submit both raw and curated data the NCBI repositories.
* Carry out a complete metagenomics data analysis, including the report of the biological findings.

## Session 9 - FAIR data and reproducibility

## Material

* **Session:** [:fontawesome-solid-file-pdf: Session slides](../assets/pdf/Intro%20MDA-Day5-Part1.pdf){: .md-button }

### Why FAIR?

Everything we did in this course produces valuable research objects: raw reads, assemblies, MAGs, taxonomy tables, figures, etc. If these outputs cannot be **found, accessed, interoperated and reused** by the research community, much of their value is lost. The **FAIR** principles provide the checklist:

* **F**indable: with a persistent identifier (DOI/accession) and rich, indexed metadata.
* **A**ccessible: deposited in a trusted repository, retrievable with standard protocols (and clear license/conditions).
* **I**nteroperable: using community formats and controlled vocabularies/ontologies, so that both humans and machines can process the data.
* **R**eusable: with clear provenance (how the data were generated and processed) and a license that states how it may be reused.

### FAIR in practice: repositories and metadata

For metagenomics, this translates into concrete actions:

* **Raw reads**: deposited in INSDC repositories ([NCBI SRA](https://www.ncbi.nlm.nih.gov/sra), [ENA](https://www.ebi.ac.uk/ena), DDBJ) at submission of the associated BioProject/BioSample. Note that the dataset used in this course ([PRJNA256993](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA256993) on [Day 2](day2_targeted_amplification_analysis.md)) is public exactly because it was submitted with proper metadata.
* **Metadata**:follow the community checklists: **MIxS** (Genomic Standards Consortium) packages for environment, host and sampling metadata; **MIMAG** for MAG reports (quality metrics). Metadata richness is what makes the difference between a reusable dataset and an uninterpretable FASTQ file.
* **Assemblies and MAGs**: deposited in GenBank/MGnify with their quality statistics and taxonomy.
* **Analysis outputs and code**: supplementary tables with DOIs (e.g. Zenodo, connected to a GitHub release), plus the exact commands and parameters used.

### Reproducible analyses: workflows and containers

FAIR data alone is not enough: the **process** that turns reads into results must also be reproducible. Throughout this course we have (implicitly and explicitly) relied on the two pillars of computational reproducibility in bioinformatics:

* **Workflow managers** — tools such as [Nextflow](https://www.nextflow.io/) (used within TaxoFlow on [Day 3](day3_reference_based_compositional_analysis.md) and for nf-core/mag on [Day 4](day4_mag_annotation.md)) capture the entire analysis, including every tool, parameter, input and output in a single executable, versioned definition that can run on a laptop, a cluster or the cloud. Communities such as [nf-core](https://nf-co.re/) provide curated, peer-reviewed workflows that you can reuse.
* **Containers** — Docker/Singularity images freeze the software environment (OS libraries, tool versions) of each step, so the same workflow produces the same result anywhere, forever. nf-core pipelines, for example, ship every process as a container.

A practical recipe for a FAIR and reproducible metagenomics analysis:

1. Version your workflow and scripts in **Git** (e.g. on GitHub), with a README that states how to run them.
2. Run them with a **workflow manager** (Nextflow, Snakemake) using **containers/conda** environments pinned to specific versions.
3. Record the **parameters and database versions** (e.g. the exact Kraken2 database and Silva release used above).
4. Deposit raw data and derived objects in **trusted repositories** with rich metadata and persistent identifiers.
5. Report according to community standards (**MIxS**, **MIMAG**), and license your code and data so others can reuse them.

## Session 10 - Real-world project

**Objective:** You are equipped now with the knowledge and skills to perform metagenomics data analysis from different perspectives. This task is designed to evaluate how you would carry out the analysis in a real-world environment.

### Data

You are provided with a small metagenomic dataset of yogurt and dietary supplement samples, accessible
via NCBI BioProject [PRJNA603575](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA603575). You will find as well the metadata for naming and grouping the samples are supplied separately in tab-delimited format (_see below_). It also contains the SRA run accessions to simplified sample names and sample types (Yogurt or Dietary_supplement), which should be used for downstream grouping and statistical comparisons.

You can download the metadata found here:

[:fontawesome-solid-file-code: Metadata](../assets/data/PRJNA603575_metadata.tsv){: .md-button }

or extract it from here:

??? abstract "Experiment metadata"
    | Run | Sample | BioSample | SampleName | SampleType |
    |---|---|---|---|---|
    | SRR11605263 | SRS6525470 | SAMN14677880 | DS001 | DietSuppl |
    | SRR11605262 | SRS6525471 | SAMN14677881 | DS002 | DietSuppl |
    | SRR11605271 | SRS6525478 | SAMN14677882 | DS003 | DietSuppl |
    | SRR11605270 | SRS6525479 | SAMN14677883 | DS004 | DietSuppl |
    | SRR11605269 | SRS6525480 | SAMN14677884 | Y001 | Yogurt |
    | SRR11605268 | SRS6525481 | SAMN14677885 | Y002 | Yogurt |
    | SRR11605267 | SRS6525482 | SAMN14677886 | Y003 | Yogurt |
    | SRR11605266 | SRS6525483 | SAMN14677887 | Y004 | Yogurt |
    | SRR11605265 | SRS6525484 | SAMN14677888 | Y005 | Yogurt |
    | SRR11605264 | SRS6525485 | SAMN14677889 | Y006 | Yogurt |
    | SRR11605261 | SRS6525472 | SAMN14677890 | Y007 | Yogurt |
    | SRR11605260 | SRS6525473 | SAMN14677891 | Y008 | Yogurt |
    | SRR11605259 | SRS6525474 | SAMN14677892 | Y009 | Yogurt |
    | SRR11605274 | SRS6525475 | SAMN14677893 | Y010 | Yogurt |
    | SRR11605273 | SRS6525476 | SAMN14677894 | Y011 | Yogurt |
    | SRR11605272 | SRS6525477 | SAMN14677895 | Y012 | Yogurt |

### Tasks

**Your task, then, is to perform the proper statistical analysis: use the output of the bioinformatics step to perform diversity analyses and differential abundance analysis between dietary supplement and yogurt samples.**

!!! info "Workflow used"
    You are receiving the output of the pipeline [TaxoFlow](https://taxoflow.work/), executed using the indexed genome `Human / GRCh38 no-alt` for Bowtie2 and the database `PlusPFP (Standard plus Refseq protozoa, fungi & plant) v04/09/2024` for Kraken2 and Bracken. The output includes the BIOM file and other intermediate files. You can find them [here](https://github.com/jeffe107/metagenomics-training/exercise).

In your analysis you should (not an exhaustive list):

1. Filters the abundance data:
    * Remove non-bacterial and non-fungal reads.
    * Check how many reads were removed.
2. Alpha diversity analyses:
    * Calculate richness and diversity indices.
    * Statistically assess differences in these indices between food types (dietary supplement vs.
yogurt).
3. Beta diversity analyses:
    * Compute Bray-Curtis dissimilarities between samples.
    * Generate a Principal Coordinate Analysis (PCoA) plot.
    * Provide statistical testing to determine whether the food types differ significantly.
4. Differential abundance analysis (DAA):
    * Identify microbial taxa significantly differing in abundance between dietary supplement and yogurt
samples, using an appropriate statistical method.
5. Biological interpretation:
    * In a short paragraph, discuss the biological interpretation of your alpha- and beta-diversity
findings. If possible include supporting references for your findings.
6. Provide the code and associated data in a manner that enables reproducibility and Open Access.

!!! warning "Use of AI"
    The use of AI is welcome and advantageous in this exercise. However, it must be used carefully, the code should be reviewed by you, and you must be prepared to defend the methodological decisions (including the choice of TaxoFlow, even though it wasn't your call).

??? tip "Import data into R"
    In case you get stuck by importing the data into R, we can help you for you to continue with the rest of the exercise:

    !!! terminal "RStudio"
        ```r title="Importing files" linenums="1"
            library(phyloseq)
            library(dplyr)

            # Importing output from TaxoFlow and metadata
            biom <- import_biom("~/Desktop/exercise/merged.biom") # Change to your local PATH
            tab_sam <- read.csv("~/Desktop/exercise/PRJNA603575_metadata.tsv", sep="\t") # Change to your local PATH

            # Process the data to generate the Phyloseq object along with the experiment information
            samples_df <- tab_sam %>% 
            tibble::column_to_rownames("Run") 
            samples = sample_data(samples_df)
            otu_inf = biom@otu_table
            colnames(otu_inf) <- gsub("-", "", colnames(otu_inf))
            phy_inf = biom@tax_table
            phy_obj <- phyloseq(otu_inf,phy_inf,samples)

            # Rename the column names in the taxonomy table
            phy_obj@tax_table@.Data <- substring(phy_obj@tax_table@.Data, 4)
            colnames(phy_obj@tax_table@.Data)<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

            tax_table_df <- as.data.frame(tax_table(phy_obj))

            # Modify the species column to include genus
            tax_table_df$Species <- paste(tax_table_df$Genus, tax_table_df$Species, sep = " ")

            # Update the taxonomy table in the phyloseq object
            tax_table(phy_obj) <- as.matrix(tax_table_df)
            phy_obj
        ```