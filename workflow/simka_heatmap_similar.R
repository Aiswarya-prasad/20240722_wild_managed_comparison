# Load required libraries

if(!require(tidyverse)){
  install.packages(pkgs = 'tidyverse', repos = 'https://stat.ethz.ch/CRAN/')
  library(tidyverse)
}

if(!require(gplots)){
  install.packages(pkgs = 'gplots', repos = 'https://stat.ethz.ch/CRAN/')
  library(gplots)
}

if(!require(RColorBrewer)){
  install.packages(pkgs = 'RColorBrewer', repos = 'https://stat.ethz.ch/CRAN/')
  library(RColorBrewer)
}

# Collect arguments

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop("Missing arguments", call.=FALSE)
} else {
  input <- args[1]
  output.heatmap <- args[2]
  output.table <- args[3] 
  dist <- args[4]
}