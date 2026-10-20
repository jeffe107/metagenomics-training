# Metagenomics Data Analysis Workshop

![Static Badge](https://img.shields.io/badge/release%20date-tbd-%230DC09D)

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

* Website is hosted at [https://sib-swiss.github.io/metagenomics-training/](https://sib-swiss.github.io/metagenomics-training/)

* Please refer to the repository [issues](https://github.com/sib-swiss/metagenomics-training/issues) for improvements/bugs for course material or the website

* Any contribution to this course material is highly appreciated :+1:. Please have a look at the [CONTRIBUTING.md](CONTRIBUTING.md) file to learn more on how to contribute

## Authors

* [Jeferyd Yepes-García](https://jeferydyepes.com/)
* Laurent Falquet

## Course content

* **Day 1 - Quality control of raw data**: quality control overview, read types, a FastQC hands-on exercise with short (Illumina), long (PacBio) and metagenomics datasets.
* **Day 2 - Targeted sequencing analysis with DADA2**: processing of targeted 16S amplicon (eDNA) reads with DADA2 (filtering, error correction, merging, chimera removal, taxonomy assignment), followed by a downstream Phyloseq analysis (rarefaction, α- and β-diversity, barplots, geographic mapping and LEfSe).
* **Day 3 - Reference-based compositional analysis**: a Nextflow-based Kraken2/Bracken pipeline, Krona plots, Pavian, and a downstream Phyloseq analysis (α- and β-diversity, ordination and networks).
* **Day 4 - MAG analysis**: integral analysis of Metagenome-Assembled Genomes (MAGs) - quality assessment with BIgMAG, metabolic potential with KEGG Decoder, functional annotation (COG, CAZymes, biosynthetic gene clusters), genome visualization with Proksee and pangenome analysis.

More material may be appended as new days are added to the course.

## How to host the website locally?

Once you have cloned the repo, you can host it on your local browser. The website is generated with [MkDocs](https://www.mkdocs.org/), with the [Material](https://squidfunk.github.io/mkdocs-material/) theme.

* Clone the repository:
	```bash
	git clone https://github.com/sib-swiss/metagenomics-training
	```

* Install MkDocs:
	```bash
	pip install mkdocs
	```

* And Material:
	```bash
	pip install mkdocs-material
	```

* Make sure you are in the repository directory and type:
	```bash
	mkdocs serve
	```

The website will be hosted on your local browser at [http://localhost:8000/](http://localhost:8000/).

## How to generate a github page?

For an automatically generated github page, you can run:

```sh
mkdocs gh-deploy
```