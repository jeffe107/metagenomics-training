# Day 4 - MAG analysis

## Learning outcomes

**After having completed this chapter you will be able to:**

* Explain how MAGs are recovered, polished and pre-processed with nf-core/mag.
* Assess the overall quality, contamination and taxonomy of MAGs with the BIgMAG dashboard (BUSCO and CheckM2).
* Infer the metabolic potential of MAGs from KEGG orthology terms with KEGG Decoder.
* Annotate and compare functional categories with COGclassifier, dbCAN3 (CAZymes) and antiSMASH (biosynthetic gene clusters).
* Visualize and annotate prokaryotic genomes with Proksee.
* Compare the gene repertoire of related MAGs with a pangenome analysis (Roary/Phandango).

## Objective

<p style='text-align: justify;'>Perform the integral analysis of Metagenome-Assembled Genomes (MAGs) recovered from the public dataset used in [Day 3](day3_reference_based_compositional_analysis.md). From the quality assessment and taxonomic annotation of the bins recovered by the assembly pipeline, to the exploration of their metabolic potential, functional annotation and pangenome, this session walks you through an end-to-end MAG analysis pipeline.</p>

## nf-core/mag

<p style='text-align: justify;'>The pipeline used to build the MAGs was [nf-core/mag](https://nf-co.re/mag/5.4.2), which integrates several tools to assemble the sequences and polish the recovered MAGs. The usage of the pipeline can be checked [here](https://nf-co.re/mag/5.4.2), and below is the workflow it follows:</p>

<figure markdown align="center">
  ![nf-core/mag workflow](../assets/images/mag_workflow.png)
  <figcaption>nf-core/mag workflow used to assemble the sequencing reads and recover the MAGs.</figcaption>
</figure>

<p style='text-align: justify;'>Since many of the downstream processes demand high computational resources, we have pre-computed some of them for you (you can find the commands to execute nf-core/mag on the **Moodle page** as `Commands to run nf-core/mag and MAGFlow`). However, we will explain step-by-step what we did or we are going to provide the references to review what's happening under the hood.</p>

## MAGFlow/BIgMAG

<p style='text-align: justify;'>[MAGFlow](https://f1000research.com/articles/13-640) is a tool designed to combine several tools to measure the quality of the bins/MAGs, as well as to taxonomically annotate them. This is the workflow:</p>

<figure markdown align="center">
  ![MAGFlow workflow](../assets/images/workflow.png)
  <figcaption>MAGFlow workflow used to assess the quality and taxonomic annotation of the MAGs.</figcaption>
</figure>

<p style='text-align: justify;'>The output of this tool is a ready-to-use file that concatenates the results, and it can then be used as input for [BIgMAG](https://f1000research.com/articles/13-640). We will not execute MAGFlow today (you can find the commands to execute MAGFlow on the **Moodle page** as `Commands to run nf-core/mag and MAGFlow`), and hence you will find the `final_df.tsv` within `intermediate/magflow/` (download and extract the folder `intermediate.tar.gz` from the **Moodle page**) to display the BIgMAG dashboard. Let's create then the environment for BIgMAG by running the following commands on the GitHub Codespaces terminal:</p>

```bash title="Set up the BIgMAG environment" linenums="1"
git clone https://github.com/jeffe107/BIgMAG
conda create -n BIgMAG --file BIgMAG/requirements.txt
conda activate BIgMAG
```

!!! warning "Whenever you are asked whether to install extra packages"
    Please say yes to all!

<p style='text-align: justify;'>Now, we are ready to execute BIgMAG to perform the exploratory analysis of the MAG overall quality and annotation:</p>

```bash title="Run the BIgMAG dashboard" linenums="1"
BIgMAG/app.py -p 8050 data/magflow/final_df.tsv
```

<p style='text-align: justify;'>You will see on the terminal a link to the dashboard or the editor will offer you an option to directly open a new tab on the browser, just click on it.</p>

<p style='text-align: justify;'>Now, it is your turn to analyze the results, use these questions to guide your thoughts:</p>

!!! question "BIgMAG - overall quality"
    Overall, the MAGs recovered within your assigned samples have good quality? What would you suggest to improve the quality of the bins that do not depict enough quality?

!!! question "BIgMAG - quality metrics"
    The results displayed by BUSCO and CheckM2 are correspondent? Why do you think this is happening?

!!! question "BIgMAG - taxonomy"
    What about the taxonomic classification, what would you report in terms of comparison among samples?

!!! question "BIgMAG - clustering"
    Were the samples clustered as you expected? What caused this behavior?

!!! question "BIgMAG - interesting cases"
    Is there any unusual sample or MAG that catches your attention? What further analysis would you propose to follow in order to dig into this special sample or MAG?

!!! warning "Stop the dashboard"
    Do not forget to stop the dashboard with `Ctrl + c`, and to deactivate the environment with:

    ```bash title="Deactivate the BIgMAG environment" linenums="1"
    conda deactivate
    ```

## KEGG Decoder

<p style='text-align: justify;'>Moving forward, it is time now to identify genomic features within the bins/MAGs. [KEGG Decoder](https://github.com/bjtully/BioData/tree/master/KEGGDecoder) is a tool that interprets the metabolic potential of the MAGs (or any genome in general) by analyzing the presence of KEGG Orthology (KO) identifiers. It maps these KOs to major metabolic pathways (biogeochemically relevant) and summarizes the completeness of each pathway based on predefined modules. Being so, we can infer the metabolic capabilities of the community, and gain insights into functional differences across samples or MAGs.</p>

<p style='text-align: justify;'>[Prokka](https://github.com/tseemann/prokka) (included within nf-core/mag) has provided us with the annotation about the presence/localization of enzymes in the bins/MAGs that are involved in a wide variety of processes. Nonetheless, Prokka uses EC numbers to describe such genomic features, and as a result, we needed to transform these annotations to K numbers and merge them into a single file (this was already done for you, the required files to perform these translations can be found within the directory stored on the **Moodle page** as `additional_files.tar.gz`). Now, let's install KEGG Decoder:</p>

```bash title="Install KEGG Decoder" linenums="1"
conda create -n keggdecoder python=3.6
conda activate keggdecoder
python3 -m pip install KEGGDecoder
```

!!! warning "Whenever you are asked whether to install extra packages"
    Please say yes to all!

<p style='text-align: justify;'>Now, we are ready to execute this software to determine the presence and completeness of the annotated metabolic pathways:</p>

```bash title="Run KEGG Decoder" linenums="1"
KEGG-decoder --input data/k_numbers/megahit_k_numbers.tsv --output kegg_output.tsv --vizoption static
```

<p style='text-align: justify;'>This command will create a heatmap where you can perform a comparison across MAGs. In this case, we are analyzing the MAGs obtained with the assembler MEGAHIT. It should look like this:</p>

<figure markdown align="center">
  ![KEGG Decoder heatmap](../assets/images/test.png)
  <figcaption>KEGG Decoder heatmap showing the presence and completeness of metabolic pathways across the MEGAHIT-assembled MAGs.</figcaption>
</figure>

<p style='text-align: justify;'>On the y axis you will see just random letters per row, these are the corresponding names MEGAHIT has assigned to the contigs in each MAG. Below you can find a mapping file between the name of the bin/MAG and the contig name ([download](../assets/data/summary_contigs.tsv)):</p>

| Bin | Contig name |
|-----|-------------|
| MEGAHIT-MaxBin2-ERR2143759.001.tsv | GIFAEIHF |
| MEGAHIT-MaxBin2-ERR2143759.002.tsv | CIJAKMNO |
| MEGAHIT-MaxBin2-ERR2143759.003.tsv | IACMDEGO |
| MEGAHIT-MaxBin2-ERR2143759.004.tsv | FKGEEGMK |
| MEGAHIT-MaxBin2-ERR2143759.005.tsv | NBLBKFGP |
| MEGAHIT-MaxBin2-ERR2143760.001.tsv | CDMCPOJD |
| MEGAHIT-MaxBin2-ERR2143760.002.tsv | PJEJNNEC |
| MEGAHIT-MaxBin2-ERR2143760.003.tsv | ELBPOCGG |
| MEGAHIT-MaxBin2-ERR2143760.004.tsv | LPDNOCIE |
| MEGAHIT-MaxBin2-ERR2143771.001.tsv | KEBGIGOI |
| MEGAHIT-MaxBin2-ERR2143771.002.tsv | GNGAKHBD |
| MEGAHIT-MaxBin2-ERR2143772.001.tsv | NGILDPFJ |
| MEGAHIT-MaxBin2-ERR2143772.002.tsv | BNBNNLFK |
| MEGAHIT-MaxBin2-ERR2143773.001.tsv | OLADMAGF |
| MEGAHIT-MaxBin2-ERR2143773.002.tsv | EPNPGCJF |
| MEGAHIT-MaxBin2-ERR2143773.003.tsv | OGHEFGBH |
| MEGAHIT-MetaBAT2-ERR2143759.1.tsv | EOMIBOKB |
| MEGAHIT-MetaBAT2-ERR2143759.2.tsv | IMFCMJNE |
| MEGAHIT-MetaBAT2-ERR2143759.3.tsv | JPIOBONM |
| MEGAHIT-MetaBAT2-ERR2143759.4.tsv | IIKMGHIJ |
| MEGAHIT-MetaBAT2-ERR2143759.5.tsv | FJMNHDHH |
| MEGAHIT-MetaBAT2-ERR2143759.6.tsv | HNFOOLAK |
| MEGAHIT-MetaBAT2-ERR2143759.7.tsv | KKNACEHJ |
| MEGAHIT-MetaBAT2-ERR2143760.1.tsv | JMOHPCEH |
| MEGAHIT-MetaBAT2-ERR2143760.2.tsv | JNLKELBO |
| MEGAHIT-MetaBAT2-ERR2143760.3.tsv | BPOICHMO |
| MEGAHIT-MetaBAT2-ERR2143760.4.tsv | OJEIIDGO |
| MEGAHIT-MetaBAT2-ERR2143760.5.tsv | OINBHHLA |
| MEGAHIT-MetaBAT2-ERR2143771.1.tsv | EONPCBDF |
| MEGAHIT-MetaBAT2-ERR2143771.2.tsv | OJCBHCCG |
| MEGAHIT-MetaBAT2-ERR2143771.3.tsv | JADMGMAL |
| MEGAHIT-MetaBAT2-ERR2143771.4.tsv | KKGCJLIG |
| MEGAHIT-MetaBAT2-ERR2143772.1.tsv | KHPHFBJJ |
| MEGAHIT-MetaBAT2-ERR2143772.2.tsv | NECFOIJH |
| MEGAHIT-MetaBAT2-ERR2143772.3.tsv | KDCCMDHF |
| MEGAHIT-MetaBAT2-ERR2143772.4.tsv | FDKDGLCP |
| MEGAHIT-MetaBAT2-ERR2143773.1.tsv | IJLNKNIL |
| MEGAHIT-MetaBAT2-ERR2143773.2.tsv | MACHMHHK |
| MEGAHIT-MetaBAT2-ERR2143773.3.tsv | IGFKCGOH |
| MEGAHIT-MetaBAT2-ERR2143773.4.tsv | LLMAHDOM |

<p style='text-align: justify;'>Alternatively, you can upload the K numbers (`data/k_numbers/megahit_k_numbers.tsv` from Codespaces) to [KEGG Mapper – Reconstruct](https://www.genome.jp/kegg/mapper/reconstruct.html) to visualize the presence/absence of all KEGG pathways. Please keep in mind that this file contains all K numbers for MEGAHIT-assembled bins/MAGs, thus what you obtain is a community-level analysis. This is one example of how it should look like:</p>

<figure markdown align="center">
  ![KEGG Mapper – Reconstruct output](../assets/images/kegg_mapper.png)
  <figcaption>KEGG Mapper – Reconstruct output depicting the presence/absence of all KEGG pathways at the community level.</figcaption>
</figure>

!!! warning "Deactivate the environment"
    Do not forget to deactivate the environment with:

    ```bash title="Deactivate the keggdecoder environment" linenums="1"
    conda deactivate
    ```

## COG Annotation

<p style='text-align: justify;'>Clusters of Orthologous Genes (COGs) are groups of genes from different organisms that evolved from a common ancestral gene and retain the same function. To explore these genes, we rely on the amino-acid sequences of the coding regions provided by Prokka (_.faa_ files). The next task, then, will be to install a [COGclassifier](https://github.com/moshi4/COGclassifier) in order to detect Clusters of Ortholog Genes:</p>

```bash title="Install COGclassifier" linenums="1"
conda create -n cogclassifier -c conda-forge -c bioconda cogclassifier
conda activate cogclassifier
```

!!! warning "Whenever you are asked whether to install extra packages"
    Please say yes to all!

<p style='text-align: justify;'>This tool _automatically performs the processes from searching query sequences into the COG database, to annotation and classification of gene functions, to generation of publication-ready figures_. However, the input for this tool is only one genome each time. This does not mean that we cannot analyze all the bins/MAGs; we could simply execute the software for each of them and integrate the data afterwards. To launch the tool just run on the terminal:</p>

```bash title="Run COGclassifier" linenums="1"
COGclassifier -i data/faas/MEGAHIT-MetaBAT2-ERR2143759.7.faa -o cog_annotation --download_dir ./cog_database
```

<p style='text-align: justify;'>Once it has finished, inside `cog_annotation` you will see the output represented as tables featuring counts, summary and annotations, as well as interesting figures showcasing the proportion of the different COG categories, just like this:</p>

<figure markdown align="center">
  ![COGclassifier output](../assets/images/cog.png)
  <figcaption>COGclassifier figure showcasing the proportion of the different COG categories for one MAG.</figcaption>
</figure>

!!! question "COG categories"
    What are the categories that are more representative of this MAG? Do you see any odd results?

!!! question "COG comparisons"
    If you decide to analyze all of the MAGs, do you think it would be a fair comparison just using the raw counts obtained by the software? If not, what strategy would you propose to proceed further?

!!! warning "Deactivate the environment"
    Do not forget to deactivate the environment with:

    ```bash title="Deactivate the cogclassifier environment" linenums="1"
    conda deactivate
    ```

## dbCAN3

<p style='text-align: justify;'>CAZymes (Carbohydrate-Active enZymes) are enzymes involved in the breakdown, biosynthesis, or modification of carbohydrates and glycoconjugates. They play a crucial role in processing complex carbohydrates such as cellulose, hemicellulose, starch, and chitin, among others. To detect the presence of this kind of enzymes, we are going to use the tool [dbCAN3](https://bcb.unl.edu/dbCAN2/blast.php), which is an automated web server designed to run the software and provide the results.</p>

<p style='text-align: justify;'>To execute the tool, you just need to download the file from **Codespaces** `data/faas/MEGAHIT-MetaBAT2-ERR2143759.7.faa`, upload it to the server (click [here](https://pro.unl.edu/dbCAN2/blast.php)), submit the job and wait for the results.</p>

!!! warning "Tools included in the execution"
    You can select the tools to include in the execution (HMMER, DIAMOND, dbCAN-sub, CGCFinder). However, please keep in mind that the more tools you include, the longer the run.

<p style='text-align: justify;'>The output should look like this (a preview of the first rows is shown below, and the full table is available as a [downloadable file](../assets/data/dbcan_output.csv)):</p>

| GeneID | EC | HMMER | dbCAN_sub | DIAMOND | SignalP | NofTools |
|--------|----|-------|-----------|---------|---------|----------|
| OJEIIDGO_00076 | None | GH2(51-929) | None | None | Y(1-22) | 1 |
| OJEIIDGO_00091 | None | GH20(151-512) | None | None | Y(1-20) | 1 |
| OJEIIDGO_00143 | None | CE1(19-254) | None | None | N | 1 |
| OJEIIDGO_00156 | None | GT41(135-688) | None | None | N | 1 |
| OJEIIDGO_00303 | None | GH16_3(35-227) | None | None | N | 1 |
| OJEIIDGO_00341 | None | GH73(158-287) | None | None | N | 1 |
| OJEIIDGO_00345 | None | GH177(38-426) | None | None | Y(1-35) | 1 |
| OJEIIDGO_00582 | None | AA3(3-568) | None | None | N | 1 |
| OJEIIDGO_00690 | None | GH109(41-201) | None | None | N | 1 |
| OJEIIDGO_00879 | None | GH73(136-324) | None | None | N | 1 |

<p style='text-align: justify;'>Visit the [CAZy](https://www.cazy.org/) website if you want to know more about enzyme families and classes reported by dbCAN3. Similar to the COGclassifier, we can perform the analysis for all the MAGs, and integrate the data afterwards to achieve an overall comparison.</p>

## antiSMASH

<p style='text-align: justify;'>[antiSMASH](https://antismash.secondarymetabolites.org) (antibiotics & Secondary Metabolite Analysis SHell) is a tool that detects and analyzes biosynthetic gene clusters (BGCs) in microbial genomes. These clusters are groups of co-located genes that together encode the machinery to produce secondary metabolites—specialized compounds that are not essential for basic cellular functions.</p>

<p style='text-align: justify;'>Analogous to dbCAN3, we will use the web server to annotate one of the MAGs recovered by the pipeline. Upload the file `intermediate/gbks/MEGAHIT-MetaBAT2-ERR2143759.7.gbk` (download and extract the folder `intermediate.tar.gz` from the **Moodle page**) to [antiSMASH](https://antismash.secondarymetabolites.org) and wait for the results.</p>

<p style='text-align: justify;'>If it takes some time, we have performed this step for you and we stored the results within the same `intermediate.tar.gz` file, at `intermediate/antiSMASH/index.html`.</p>

## Proksee

<p style='text-align: justify;'>[Proksee](https://proksee.ca/) is an interactive web-based tool for visualizing, annotating, and analyzing prokaryotic genomes. Using this tool we can visualize genomes as circular or linear maps, annotate them, customize the visualization and export high-level and detailed-oriented figures.</p>

<p style='text-align: justify;'>Same as with previous applications, you just need to upload the annotated genome in FASTA or GenBank format. Here, we are going to leverage the files produced by Prokka, and hence you just need to upload the file (download and extract the folder `intermediate.tar.gz` from the **Moodle page**) `intermediate/gbks/MEGAHIT-MetaBAT2-ERR2143759.7.gbk` to [Proksee](https://proksee.ca/).</p>

<p style='text-align: justify;'>You can customize the display, add features, re-annotate the genome among many other functionalities. Unfortunately, it processes only one genome at a time. From here it's up to your goals/motivation to take advantage of this tool.</p>

<p style='text-align: justify;'>You should be seeing this example:</p>

<figure markdown align="center">
  ![Proksee circular genome map](../assets/images/proksee.png)
  <figcaption>Proksee circular map visualizing the annotated genome of one of the MAGs.</figcaption>
</figure>

## Pangenome

<p style='text-align: justify;'>For the past tools we have not analyzed just a random MAG: we selected this one given that its taxonomic annotation is shared across different samples (from the BIgMAG exercise), and therefore it is interesting to establish the similarities/differences among these MAGs to enable a pangenome analysis. In this case, we study the entire gene repertoire of related MAGs (same species or genus); for our analysis we have selected some MAGs based on GTDB classification.</p>

<p style='text-align: justify;'>Next, we run the tool Roary that determines the core genome, accessory genome and unique genome. We have pre-computed the results for you following the [tutorial](https://sanger-pathogens.github.io/Roary/) presented by the developers of the tool, and now we are going to visualize the results using [Phandango](https://jameshadfield.github.io/phandango/).</p>

<p style='text-align: justify;'>Drag and drop the files (download and extract the folder `intermediate.tar.gz` from the **Moodle page**) `intermediate/pangenome/JAAUTG01/workshop.newick` and `intermediate/pangenome/JAAUTG01/gene_presence_absence.csv` to the [Phandango web server](https://jameshadfield.github.io/phandango/) and visualize the results.</p>

<p style='text-align: justify;'>It should look like this:</p>

<figure markdown align="center">
  ![Phandango pangenome visualization](../assets/images/phandango.png)
  <figcaption>Phandango output visualizing the Roary pangenome analysis of the selected MAGs.</figcaption>
</figure>

<p style='text-align: justify;'>Usually, this analysis is carried out using a reference genome; however, given that these MAGs are not annotated at species level, and genus annotation is not informative, we do not count with a reference genome to explore the pangenome of these MAGs. You can visualize an example that includes a reference genome with the files found at `intermediate/pangenome/example_with_reference` (download and extract the folder `intermediate.tar.gz` from the **Moodle page**).</p>

## Additional resources

<p style='text-align: justify;'>If you want to explore more in detail how you can expand your analysis, please check the following resources:</p>

* [nf-core/mag pipeline documentation](https://nf-co.re/mag/)
* [MAGFlow and BIgMAG publication](https://f1000research.com/articles/13-640)
* [KEGG Decoder repository](https://github.com/bjtully/BioData/tree/master/KEGGDecoder)
* [COGclassifier repository](https://github.com/moshi4/COGclassifier)
* [dbCAN web server](https://bcb.unl.edu/dbCAN2/)
* [CAZy database](https://www.cazy.org/)
* [antiSMASH web server](https://antismash.secondarymetabolites.org)
* [Proksee web server](https://proksee.ca/)
* [Roary tutorial](https://sanger-pathogens.github.io/Roary/)
* [Phandango web server](https://jameshadfield.github.io/phandango/)