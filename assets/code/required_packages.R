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