# Day 4 - Metagenome assembly and binning

## Learning outcomes

**After having completed this chapter you will be able to:**

* Explain the principles of _de novo_ metagenomic assembly and the algorithms (overlap–layout–consensus, de Bruijn graphs) behind the common metagenome assemblers.
* Discuss the challenges that community complexity, strain microdiversity and uneven coverage pose for assembly, and how they translate into fragmented or chimeric contigs.
* Evaluate assembly quality with the standard metrics and tools (contig count, N50, MetaQUAST).
* Describe binning strategies for recovering individual genomes from assemblies, and distinguish composition-based from coverage-based (and hybrid) approaches.
* Explain how differential coverage across samples improves binning, and how assemblers and binners are combined into workflows such as nf-core/mag.
* Inspect the output of an end-to-end MAG-recovery pipeline and decode the naming scheme of the recovered bins.

## Objective

In this chapter we recover **Metagenome-Assembled Genomes (MAGs)** from the shotgun dataset introduced in [Day 3](day3_reference_based_compositional_analysis.md) (Cuatro Ciénegas, México). **Session 7** covers the theory of _de novo_ metagenome assembly; **Session 8** covers how individual genomes are carved out of the assembly (binning) and how both steps are orchestrated in a single workflow with [nf-core/mag](https://nf-co.re/mag/). The quality assessment, annotation and interpretation of the recovered MAGs is the topic of [Day 5](day5_mag_annotation_and_fair.md).

## Session 7 - Metagenomics assembly

### Principles of _de novo_ assembly

Assembly reconstructs the original genomes from millions of short (or long) reads, without a reference. In its simplest form the problem is solved by finding overlaps between reads and laying them out into longer contiguous sequences (**contigs**), a paradigm known as **overlap–layout–consensus (OLC)** — the basis of classic assemblers such as CAP3 or Canu/Flye for long reads.

Modern short-read assemblers instead build a **de Bruijn graph**: the reads are broken into k-mers (overlapping fragments of length _k_), which become nodes of a graph in which edges represent (k–1)-mer overlaps. Reads from the same genomic region then trace paths through the graph, and contigs are extracted as non-branching paths. This is the foundation of the metagenome assemblers we use in practice:

* **MEGAHIT** — compact, fast de Bruijn graph assembler designed for large and complex metagenomes (single-end or paired-end short reads). Used in the pipeline below.
* **metaSPAdes / metaViSPAdes** — multi-k de Bruijn approach with better contiguity at higher compute cost.
* **metaFlye** — Flye adapted to noisy long reads; the standard choice for Nanopore/PacBio assemblies, often combined with short reads in hybrid mode.
* **hifiasm-meta** — designed for high-accuracy long reads (PacBio HiFi).

### Challenges of metagenome assembly

Assembling a *community* is fundamentally harder than assembling a single genome:

* **Uneven coverage** — abundant organisms assemble into long contigs, while rare members remain fragmented or unassembled altogether.
* **Strain microdiversity** — closely related strains share most of their sequence; their strain-specific variants create branches in the graph and shatter assemblies (or worse, produce chimeric contigs that mix strains).
* **Community complexity** — soil or marine communities contain thousands of genomes; the more complex the community, the shorter and more tangled the contigs.
* **Uneven genome characteristics** — GC extremes and repeats further break contiguity.

The consequences are chimeric contigs, fragmented genomes and mis-binned sequences downstream — always keep this in mind when interpreting the length and quality distributions of what we will recover in Session 8.

### Assembly quality metrics

Assemblies are typically summarized with:

* **Number of contigs** and total assembly size — how fragmented is the reconstruction?
* **N50 / L50** — the contig length at which 50% of the assembly is contained in contigs at least that long (a contiguity summary; not a correctness measure!).
* **Misassemblies and correctly placed bases** — assessed by aligning the assembly to a reference with tools such as **QUAST**/**MetaQUAST** (when a close reference exists) or by mapping the reads back to the assembly to check consistency (e.g. with **Bandage**, or with the mapping statistics produced by the pipeline).
* For the recovered genomes later on, quality is judged with completeness/contamination estimates (CheckM2, BUSCO) — the topic of [Day 5](day5_mag_annotation_and_fair.md).

## Session 8 - Draft genomes (bins) and Metagenome-Assembled Genomes (MAGs)

### From contigs to genomes: binning strategies

A metagenome assembly is a soup of contigs from many organisms. **Binning** groups contigs into draft genomes — the **MAGs**. Two complementary signals drive every binner:

* **Composition** — each genome has a characteristic tetranucleotide (or longer k-mer) frequency and GC content, imprinted on all of its contigs. Composition alone struggles with short contigs and with closely related strains.
* **Coverage** — the abundance of a genome's contigs should be similar. Crucially, using **multiple samples** (differential coverage) gives each genome a coverage profile across samples that acts as a powerful fingerprint, since different organisms fluctuate independently between samples.
* **Hybrid binners** — combine both signals. The workhorses in current pipelines are **MetaBAT2** and **MaxBin2** (both run by nf-core/mag, which is why you will see two series of bins per sample below), with newer alternatives such as **VAMB**, **CONCOCT** and ensemble refinements such as **DAS Tool** that reconcile the outputs of several binners.

!!! tip "How to get better bins"
    Sequence several samples that share members of the community (differential coverage!), aim for sufficient depth ([Day 3](day3_reference_based_compositional_analysis.md)), and refine/consensus-bin the outputs of multiple binners rather than trusting a single one.

### nf-core/mag: assembly and binning in one workflow

The pipeline used to build the MAGs was [nf-core/mag](https://nf-co.re/mag/5.4.2), a Nextflow workflow that chains quality control and host-read removal, assembly (MEGAHIT and/or metaSPAdes), binning (MetaBAT2, MaxBin2), bin refinement, and downstream genome annotation and quality control (Prokka, CheckM). The usage of the pipeline can be checked [here](https://nf-co.re/mag/5.4.2), and below is the workflow it follows:

<figure markdown align="center">
  ![nf-core/mag workflow](../assets/images/mag_workflow.png)
  <figcaption>nf-core/mag workflow used to assemble the sequencing reads and recover the MAGs.</figcaption>
</figure>

Since running the pipeline demands high computational resources, we have **pre-computed it for you** on the same nine samples used in [Day 3](day3_reference_based_compositional_analysis.md) (you can find the commands to execute nf-core/mag on the **Moodle page** as `Commands to run nf-core/mag and MAGFlow`). In the hands-on below we explore its outputs, and we decode what the pipeline produced.

### Hands-on: exploring the recovered MAGs

If not done already, login to your GitHub account and open the Codespaces training environment <span style="color:red">opening this link in a new tab:</span>

[https://codespaces.new/nextflow-io/training?quickstart=1&ref=2.1.5](https://codespaces.new/nextflow-io/training?quickstart=1&ref=2.1.5){: .md-button }

The pre-computed results are placed within the `metagenomics` directory of the workspace (`data/faas/` contains the predicted proteins per MAG; download and extract the folder `intermediate.tar.gz` from the **Moodle page** for the GenBank annotations and the remaining outputs used on [Day 5](day5_mag_annotation_and_fair.md)).

#### 1. Look at the recovered MAGs

The protein FASTA files (`.faa`, produced by Prokka) are named after the sample, the assembler and the binning tool. List them and count how many MAGs were recovered per binner:

```bash title="List the MAGs recovered by the pipeline" linenums="1"
cd data/faas

# how many MAGs in total?
ls *.faa | wc -l

# how many per binning tool?
ls *.faa | grep MetaBAT2 | wc -l
ls *.faa | grep MaxBin2 | wc -l

# which samples contributed MAGs?
ls *.faa | sed 's/\..*//' | sort -u
```

!!! question "Question"
    Do both binners recover the same number of MAGs per sample? What does this tell you about the (dis)agreement between composition- and coverage-based binning?

#### 2. Decode the naming scheme

```bash title="Understand the file naming" linenums="1"
# the pattern is: ASSEMBLER-BINNER-SAMPLE.BIN.faa
# e.g. MEGAHIT-MetaBAT2-ERR2143759.7.faa
#   -> assembled with MEGAHIT, binned with MetaBAT2,
#      from sample ERR2143759, bin number 7

ls MEGAHIT-MetaBAT2-*.faa
ls MEGAHIT-MaxBin2-*.faa
```

#### 3. Inspect the content of one MAG

```bash title="Peek into a MAG" linenums="1"
# how many proteins were predicted in this MAG?
grep -c '^>' MEGAHIT-MetaBAT2-ERR2143759.7.faa

# look at the predicted proteins (amino-acid sequences)
head -20 MEGAHIT-MetaBAT2-ERR2143759.7.faa
```

Each of these files is a **draft genome**: a set of genes recovered from the community without ever culturing its organism. The GenBank version of the same MAG (with genes annotated on the contigs) is found at `intermediate/gbks/MEGAHIT-MetaBAT2-ERR2143759.7.gbk` — we will use it on [Day 5](day5_mag_annotation_and_fair.md) for visualization and secondary-metabolite detection.

#### 4. Relate bins to their contigs

MEGAHIT names contigs with short random-looking codes. The mapping between each MAG and its contig name is provided in the following file ([download](../assets/data/summary_contigs.tsv)):

```bash title="Map bins to contigs" linenums="1"
head ../summary_contigs.tsv
```

| Bin | Contig name |
|-----|-------------|
| MEGAHIT-MaxBin2-ERR2143759.001.tsv | GIFAEIHF |
| MEGAHIT-MaxBin2-ERR2143759.002.tsv | CIJAKMNO |
| MEGAHIT-MetaBAT2-ERR2143759.7.tsv | KKNACEHJ |

(A subset of the full table is shown above; the complete mapping is in the downloadable file.) This kind of bookkeeping matters when you go back from a downstream result (e.g. "this MAG degrades cellulose") to the actual sequences that support it.

### Discussion

!!! question "Questions"
    * Which samples yielded MAGs, and which did not? Given what you learned in Session 7, what could explain this?
    * Why are there two series of bins (MetaBAT2 and MaxBin2) for the same assembly? Would you trust one binner over the other, and how could you decide?
    * What are the limitations of the MAGs we just recovered (completeness, contamination, strain chimerism)? How would you quantify them?

These last questions lead us directly to [Day 5](day5_mag_annotation_and_fair.md), where we assess the quality of these MAGs, annotate their functions and discuss how to make such analyses FAIR and reproducible.

## Additional resources

If you want to explore more in detail how you can expand your analysis, please check the following resources:

* [nf-core/mag pipeline documentation](https://nf-co.re/mag/)
* [MEGAHIT repository](https://github.com/voutcn/megahit)
* [metaSPAdes documentation](https://cab.spbu.ru/software/spades/)
* [metaFlye repository](https://github.com/fenderglass/Flye)
* [MetaBAT2 publication](https://doi.org/10.1186/s13059-017-1141-7)
* [MaxBin2 publication](https://doi.org/10.1093/bioinformatics/btu031)
* [DAS Tool publication](https://doi.org/10.1186/s13059-018-1540-z)
* [QUAST / MetaQUAST](https://quast.sourceforge.net/)
* [MIMAG standards for reporting MAGs](https://doi.org/10.1038/nbt.3893)
