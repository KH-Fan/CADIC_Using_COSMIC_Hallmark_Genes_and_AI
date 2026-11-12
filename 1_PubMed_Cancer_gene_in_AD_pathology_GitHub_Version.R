###-----------------------------------------------------------------------------###
###
### This script is utilizing R packages and iterate gene names, and
### also focus on the synonyms of each gene symboles 
### by using OR Boolean parameter to generate standard PubMed search terms 
### and then retrive publications doi for downstream analysis
###
###-----------------------------------------------------------------------------###

# Please set your own working directory 
#####################################
# !!! (action needed down here) !!! #
#####################################
# setwd("working directory such as ~/Desktop/")

# Optional: keep string value as a characters 
# options(stringAsFactor=FALSE)

# loading necessary packages
library(tidyverse)
library(tidylog)
library(easyPubMed)

# read in the COSMIC Hallmark Cancer gene list from the file
# this file should be under your working directory
cosmic_genes <- read.csv("Hallmark_Cancer_Genes.txt", header=FALSE)

###-----------------------------------------------------------------------------###

# build the PubMed search term
# Alzheimer's Disease as a MeSH term
# Restrict to MeSH Major Topic.
# subheading focus on genetics and pathology
# ("Alzheimer Disease/genetics"[Majr]) AND ("Alzheimer Disease/pathology"[Majr])
# then the general search term will be
general_term <- '("Alzheimer Disease/genetics"[Majr]) AND ("Alzheimer Disease/pathology"[Majr])'

# Add a personal api key
# The api key can be obtained via MYNCBI 
# the script can be run without the api key, 
# however, it will limited the search power in terms of time and efficiency
# please replace your personal api key here if you have one
#####################################
# !!! (action needed down here) !!! #
#####################################
# apikey <- "paste your own api key"

# create an empty object to store the publications
refs <- NULL

# iterate all the gene names in the search query and retrieve the publications
for (i in 1:nrow(cosmic_genes_with_synonyms)) {
  Sys.sleep(1)
  # the current gene
  cur_gene <- cosmic_genes_with_synonyms$Gene.Symbol[i]
  
  # synonyms gene names
  # set the search strategy for the gene symbol in the title and abstract only: [TIAB]
  syno_gene_search_term <- paste0(" AND ( ",cur_gene, " [TIAB] OR ", 
                                  gsub(pattern = ",", replacement = " [TIAB] OR ", 
                                  cosmic_genes_with_synonyms$Synonyms[i]), " [TIAB] )")
  
  # the current search query 
  search_query <- paste0(general_term,  syno_gene_search_term) 
  
  # print the current gene names 
  print(cur_gene)
  
  # check if there is any hit from the search query
  ids <- epm_query(query_string = search_query, api_key = apikey)
  
  # if there are matches, extract the publication information needed for the downstream analysis
  if (as.integer(ids@meta$exp_count) > 0) {
    
    # empty the objects inside the if loop
    df <- NULL
    df_subset <- NULL
    
    # get information of the publication(s)
    papers <- epm_fetch(ids, encoding = "ASCII")
    infors <- epm_parse(papers, include_abstract = FALSE, compact_output = TRUE)
    
    # Get the processed data frame out of the object
    df <- get_epm_data(infors)
    
    # Keep only the columns you need
    slctd_fields <- c('pmid', 'doi', 'title', 'journal', 'year', 'authors')
    df_subset <- subset(df, select=slctd_fields)
    df_subset$gene <- cur_gene
    refs <- rbind(refs, df_subset)
  }
}

# Today's date
todayday <- format(Sys.Date(), "%Y_%m_%d")

#####################################
# !!! (action needed down here) !!! #
#####################################
# Output directory
# Please set your own output directory if it's different than working directory
# out_dir <- "~/Desktop"

# Assign the file name with the information of today's day 
outfl <- paste0(out_dir,"/COSMIC_Genes_in_AD_Publications_with_gene_and_synonyms_in_Title_Abstract_",todayday,".txt")

# Write the outcome into a tab deliminator text file 
write.table(refs, outfl, quote=FALSE, row.names = FALSE, sep = "\t")

###-----------------------------------------------------------------------------###
### Summary Statistics
###-----------------------------------------------------------------------------###

# Read in file if generated previously
# refs <- read.csv(outfl, header=TRUE, sep="\t")

# numbers of publications in total
nrow(refs)
# 1038

# numbers of genes with publication in AD pathology publication
length(unique(refs$gene))
# 155

# Distributions of the publication for each gene
table(table(refs$gene))
# 116 genes has only one AD pathology paper, and the BRCA2 and FANCD2 gene has 116 AD pathology papers each

#####################################
# !!! (action needed down here) !!! #
#####################################
# Save the results
# Please rename your output file
# out_fl_name <- "gene_vs_publicaiton.txt"
x <- as.data.frame(table(refs$gene))
sum_pub_count_by_gene <- x[order(x$Freq),]
names(sum_pub_count_by_gene) <- c("gene", "freq")
write.table(sum_pub_count_by_gene, out_fl_name, sep="\t",
            row.names = FALSE, quote = FALSE)

###-----------------------------------------------------------------------------###
### Figure 1 in the manuscript
###-----------------------------------------------------------------------------###

# Make a table for the figure legend
sum_pub_gene_res <- sum_pub_count_by_gene %>%
  arrange(gene) %>%                     # alphabetize within each group
  group_by(freq) %>%
  summarise(genes = paste0(gene, collapse = "; ")) %>%
  arrange(freq)

# Separate each groups when the gene list is long
g1 <- sum_pub_gene_res$genes[1]
g2 <- sum_pub_gene_res$genes[2]
g3 <- sum_pub_gene_res$genes[3]
g4 <- sum_pub_gene_res$genes[4]

# Wrap to a fixed character width (e.g., 60 characters per line)
wrapped_g1 <- strwrap(g1, width = 100)
wrapped_g2 <- strwrap(g2, width = 100)
wrapped_g3 <- strwrap(g3, width = 100)
wrapped_g4 <- strwrap(g4, width = 100)

# Collapse into a single string with line breaks
wrapped_text_g1 <- paste(wrapped_g1, collapse = "\n")
wrapped_text_g2 <- paste(wrapped_g2, collapse = "\n")
wrapped_text_g3 <- paste(wrapped_g3, collapse = "\n")
wrapped_text_g4 <- paste(wrapped_g4, collapse = "\n")

# Visualize such distribution
# find the ranges of the data
bk <- range(x$Freq)
# set the range of plotting 
mi <- bk[1]-0.5
ma <- bk[2]+0.5

#####################################
# !!! (action needed down here) !!! #
#####################################
# Save the results
# Please name your output file
# pdf_name <- "Figure_1.pdf"

# plot the histogram
pdf(pdf_name, height = 9, width = 12)
par(mai=c(0,0,0,0), mar=c(4,4,2,1))
hist(table(add_refs$gene), breaks = seq(mi, ma, 1), 
     las=1, xlab="Numbers of Publications",
     ylab="Gene Counts", main="Numbers of Alzheimer's Disease Publications per Gene",
     axes = FALSE, ylim=c(0, 75))
axis(1, at=seq(0,120,5), cex.axis=0.8)
axis(2, at=seq(0,70,5), cex.axis=0.8, las=1)

# annotate the gene symbol with large number of publications

text(x=110, y=3, labels="BRCA2; FANCD2 (116)", cex=0.8, font=3)
text(x=84, y=3, labels="SET (81)", cex=0.8, font=3)
text(x=59, y=3, labels="ROS1 (58)", cex=0.8, font=3)
text(x=49, y=3, labels="AKT1 (46)", cex=0.8, font=3)
text(x=43, y=7, labels="RASA1 (43)", cex=0.8, font=3)
segments(x0=43, y0=6.4, x1=43, y1=1.2)
text(x=38, y=5, labels="PER1 (37)", cex=0.8, font=3)
segments(x0=36, y0=4.4, x1=37, y1=1.2)
text(x=32, y=3, labels="ELL (32)", cex=0.8, font=3)
text(x=25, y=5, labels="MTOR; PIK3CA (28)", cex=0.8, font=3)
segments(x0=25, y0=4.4, x1=28, y1=2.2)

# List the genes grouped by numbers of the publications
mtext(wrapped_text_g1, side=4, line=-37, las=1, cex=0.8, font=3, at=64)
mtext(wrapped_text_g2, side=4, line=-37, las=1, cex=0.8, font=3, at=57)
mtext(wrapped_text_g3, side=4, line=-37, las=1, cex=0.8, font=3, at=53)
mtext(wrapped_text_g4, side=4, line=-37, las=1, cex=0.8, font=3, at=50)
mtext(sum_pub_gene_res$genes[5], side=4, line=-37, las=1, cex=0.8, font=3, at=47)
mtext(sum_pub_gene_res$genes[6], side=4, line=-37, las=1, cex=0.8, font=3, at=44)
mtext(sum_pub_gene_res$genes[7], side=4, line=-37, las=1, cex=0.8, font=3, at=41)
mtext(sum_pub_gene_res$genes[8], side=4, line=-37, las=1, cex=0.8, font=3, at=38)
mtext(sum_pub_gene_res$genes[9], side=4, line=-37, las=1, cex=0.8, font=3, at=35)
mtext(sum_pub_gene_res$genes[10], side=4, line=-37, las=1, cex=0.8, font=3, at=32)
mtext(sum_pub_gene_res$genes[11], side=4, line=-37, las=1, cex=0.8, font=3, at=29)
mtext(sum_pub_gene_res$genes[12], side=4, line=-37, las=1, cex=0.8, font=3, at=26)
mtext(sum_pub_gene_res$genes[13], side=4, line=-37, las=1, cex=0.8, font=3, at=23)
mtext(sum_pub_gene_res$genes[14], side=4, line=-37, las=1, cex=0.8, font=3, at=20)
mtext(sum_pub_gene_res$genes[15], side=4, line=-37, las=1, cex=0.8, font=3, at=17)

# annotated the number of publication(s)
mtext("1 Publication", side=4, line=-47, las=1, cex=0.8, at=64)
mtext("2 Publications", side=4, line=-47, las=1, cex=0.8, at=57)
mtext("3 Publications", side=4, line=-47, las=1, cex=0.8, at=53)
mtext("4 Publications", side=4, line=-47, las=1, cex=0.8, at=50)
mtext("5 Publications", side=4, line=-47, las=1, cex=0.8, at=47)
mtext("6 Publications", side=4, line=-47, las=1, cex=0.8, at=44)
mtext("7 Publications", side=4, line=-47, las=1, cex=0.8, at=41)
mtext("8 Publications", side=4, line=-47, las=1, cex=0.8, at=38)
mtext("9 Publications", side=4, line=-47, las=1, cex=0.8, at=35)
mtext("12 Publications", side=4, line=-47, las=1, cex=0.8, at=32)
mtext("13 Publications", side=4, line=-47, las=1, cex=0.8, at=29)
mtext("15 Publications", side=4, line=-47, las=1, cex=0.8, at=26)
mtext("16 Publications", side=4, line=-47, las=1, cex=0.8, at=23)
mtext("17 Publications", side=4, line=-47, las=1, cex=0.8, at=20)
mtext("19 Publications", side=4, line=-47, las=1, cex=0.8, at=17)

# annotated the number of genes
mtext("[68]", side=4, line=-40, las=1, cex=0.8, at=64)
mtext("[29]", side=4, line=-40, las=1, cex=0.8, at=57)
mtext("[13]", side=4, line=-40, las=1, cex=0.8, at=53)
mtext("[8]", side=4, line=-40, las=1, cex=0.8, at=50)
mtext("[4]", side=4, line=-40, las=1, cex=0.8, at=47)
mtext("[5]", side=4, line=-40, las=1, cex=0.8, at=44)
mtext("[4]", side=4, line=-40, las=1, cex=0.8, at=41)
mtext("[3]", side=4, line=-40, las=1, cex=0.8, at=38)
mtext("[2]", side=4, line=-40, las=1, cex=0.8, at=35)
mtext("[2]", side=4, line=-40, las=1, cex=0.8, at=32)
mtext("[1]", side=4, line=-40, las=1, cex=0.8, at=29)
mtext("[2]", side=4, line=-40, las=1, cex=0.8, at=26)
mtext("[1]", side=4, line=-40, las=1, cex=0.8, at=23)
mtext("[2]", side=4, line=-40, las=1, cex=0.8, at=20)
mtext("[1]", side=4, line=-40, las=1, cex=0.8, at=17)

dev.off()

###-----------------------------------------------------------------------------###
### End of the Script
###-----------------------------------------------------------------------------###

