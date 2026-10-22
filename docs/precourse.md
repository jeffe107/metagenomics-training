## Background knowledge

As is stated in the course prerequisites, we expect participants to have a basic understanding of working with the command line on UNIX-based systems and with RStudio, a [GitHub](https://github.com/signup) account and a [GitHub Codespaces](https://github.com/features/codespaces) account.

* **UNIX skills:** test them with the quiz [here](https://docs.google.com/forms/d/e/1FAIpQLSd2BEWeOKLbIRGBT_aDEGPce1FOaVYBbhBiaqcaHoBKNB27MQ/viewform?usp=sf_link). If you don't have experience with the UNIX command line, follow this [online UNIX tutorial](https://edu.sib.swiss/pluginfile.php/2878/mod_resource/content/4/couselab-html/content.html).

* **R skills:** be comfortable with RStudio, how to install a library, matrix and data frame manipulation, and importing/exporting data from text files. If in doubt, evaluate your R skills with [this quiz](https://www.coursera.org/resources/r-programming-skill-assessment).

### R packages

You can pre-install the R/Bioconductor packages used during the course (Phyloseq ecosystem, DADA2, vegan, microbiomeMarker, etc.) with the script below. This is optional: the packages are also introduced and installed during the workshop.

[:fontawesome-solid-file-code: Download the required packages script](assets/code/required_packages.R){: .md-button }

Or, copy and paste:

!!! terminal "RStudio"
    ```R title="Required packages" linenums="1"
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

    # These might require some recompilation and ask for additional packages, please answer Yes to the questions. 

    # You can then check for each library if they load with these commands: 

    library(dada2, quietly = TRUE) 
    library(phyloseq, quietly = TRUE) 
    library(Biostrings, quietly = TRUE) 
    library(ggplot2, quietly = TRUE) 
    library(magrittr, quietly = TRUE) 
    library(vegan, quietly = TRUE) 
    library(gridExtra, quietly = TRUE) 
    library(phylogeo, quietly = TRUE) 


    #problem if microbiomeMarker does not exist for Bioc 3.23... 

    install.packages("remotes") 
    remotes::install_github("yiluheihei/microbiomeMarker") 

    library(microbiomeMarker, quietly = TRUE) 
    ```

## Software

### OS

This is an OS-agnostic course that requires only a laptop, a modern browser and a GitHub Codespaces account.

### Code editor

You will be provided with a link to create the workspace in GitHub Codespaces automatically. This link will open VS Code in your browser, so we expect you to be familiar with the layout and basic functionalities of VS Code. Otherwise, please check this quick [tutorial](https://code.visualstudio.com/docs/getstarted/getting-started) before the course to understand where everything is. Either way, there will be a short introduction to VS Code at the beginning of the online session.

This is the link to GitHub Codespaces:

[![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/sib-swiss/metagenomics-training?quickstart=1)

!!! warning "Create your Codespaces account"
    We use the free tier of [GitHub Codespaces](https://github.com/features/codespaces). Please create your account before the course, but **do not start using the Codespace yet** — the free tier has limited computing hours and storage, and they can run out quickly.

#### Video tutorial

You can find a (cool) video tutorial to learn about VS Code:

<iframe width="560" height="315" src="https://www.youtube.com/embed/1ZfO149BJvg" title="" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture; web-share" allowFullScreen><a href="https://www.ivatech.dev" style="display:none;">website development</a></iframe>

### Back-up plan

If you exceed the Codespace resources, we have prepared a solution using [CodeSandbox](https://codesandbox.io/). These are the steps to follow:

1. Go to the website: [CodeSandbox](https://codesandbox.io/).
2. Sign in using your GitHub or your Google Account (top right corner).
3. Once signed-in, click on `+ Create` (top right corner).
4. Select `Docker CodeSandbox`.
5. On the `Configure window` leave everything as default and click on Create Devbox.
6. Wait until the microVM is created.
7. Now it looks like VS Code. Create a new terminal then.
8. You are ready to run the commands from the course material.