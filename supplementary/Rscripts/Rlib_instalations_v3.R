#! /usr/bin/env Rscript
install.packages("easypackages")

# for magicians
install.packages(c("devtools", "docopts", "RColorBrewer","htmltools", "remotes"))

#transform data
install.packages(c("reshape2", "broom", "tidyverse","scales",  "data.table"))# "grid" is depreciated
#c("stringr", "readr", "dplyr", "tidyr", "tibble", "ggplot2") # part of tidyverse

#biological packages
install.packages(c("fossil","vegan", "survival", "ape", "limma"))
remotes::install_github("MadsAlbertsen/ampvis2", force=T)
## info fo fossil/chao1: https://www.rdocumentation.org/packages/fossil/versions/0.4.0/topics/chao1

# tools for Bioc-magicians (old/new version packages coalitions, expect problems, its a mess)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager") # there might be problems with newer versions and older are not supported
BiocManager::install("biomformat")
BiocManager::install("S4Vectors")






