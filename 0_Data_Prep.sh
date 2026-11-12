#!/bin/bash

###-----------------------------------------------------------------------------###
### This script is to generate files for gathering publications
### Also for the input files for AI analysis
###-----------------------------------------------------------------------------###

# Set up the working directory where the data is downloaded from COSMIC
# For example: 
# cd ~/Desktop

###
### Untar the COSMIC Hallmark Genes
###

# Current version 104, released on 19th May 2026, which was downloaded in Aug 2026
tar -vxf Cosmic_CancerGeneCensusHallmarksOfCancer_Tsv_v104_GRCh38.tar 

# Generate the list of genes and their synonymous gene IDs for R scripts to simultaneously get publications for each cancer hallmark gene
zless Cosmic_CancerGeneCensusHallmarksOfCancer_v104_GRCh38.tsv.gz | awk '{print $1}' | sort | uniq > Hallmark_Cancer_Genes.txt
# Get the list of PubMed IDs for each gene and its description of the relationship between gene and cancer
zless Cosmic_CancerGeneCensusHallmarksOfCancer_v104_GRCh38.tsv.gz | awk '{print $1, $4, $7}'  Hallmark_Cancer_Genes_and_Refs.txt

