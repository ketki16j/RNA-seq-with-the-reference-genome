# Install biomaRt package
source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

# Load biomaRt library
library("biomaRt")

# Convert final results .csv file into .txt file
results_csv <- "Gmax_DESeq2-replaceoutliers-results.csv"
write.table(read.csv(results_csv), gsub(".csv",".txt",results_csv))
results_txt <- "Gmax_DESeq2-replaceoutliers-results.txt"

#biomaRt is for http only so to get around this you need to add this line
options(RCurlOptions=list(followlocation=TRUE, postredir=2L))
# List available datasets in database
# Then select dataset to use
ensembl = useMart("phytozome_mart", host="phytozome.jgi.doe.gov",path="/biomart/martservice", dataset="phytozome")
listDatasets(ensembl)
# Check the database for entries that match the IDs of the differentially expressed genes from the results file
a <- read.table(results_txt, head=TRUE)

b <- getBM(attributes=c('transcript_id',
 'chr_name1',
 'gene_name1',
 'gene_chrom_start',
 'gene_chrom_end',
 'kog_desc'),
 filters=c('pac_transcript_id'),
 values=a$X,
 mart=ensembl)
biomart_results <- "Gmax_DESeq2-Biomart.txt"
write.table(b,file=biomart_results)
