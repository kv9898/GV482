rm(list=ls())

# setup
need <- c('tidyverse','modelsummary','haven','fixest','kableExtra',"TAM","car") # list packages needed
have <- need %in% rownames(installed.packages()) # checks packages you have
if(any(!have)) install.packages(need[!have]) # install missing packages
invisible(lapply(need, library, character.only=T))

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data <- read_dta("PopulismTerrorism.dta")
data <- as_factor(data)

sink("output/23850_log.txt")
1:10
sink()
