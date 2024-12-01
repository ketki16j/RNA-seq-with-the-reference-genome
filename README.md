# RNA-seq-with-the-reference-genome

This repository is a tutorial for differential expression analysis using RNA-Seq data
We will use a subset of gene expression data from Reid et al. 2016, a study of evolutionary adaptation to industrial pollutants in the coastal fish Fundulus heteroclitus. The study consisted of an experimental exposure of embryonic fish from eight populations to the toxicant PCB-126. Here we will use exposed and control samples from two populations.
The samples we will be using are described by the following accession numbers; SRR391535, SRR391536, SRR391537, SRR391538, SRR391539, and SRR391541.

We will be going through quality control of the reads, alignment of the reads to the reference genome, conversion of the files to raw counts, analysis of the counts with DeSeq2, and finally annotation of the reads using Biomart.

```1. The script for downloading .SRA files and converting them to fastq can be found in ~/scripts/fastq.sh```
```2. Quality Control on the Reads Using Sickle```: Step one is to perform quality control on the reads using Sickle. We are using unpaired reads, as indicated by the “se” flag in the script quality_control.sh. The -f flag designates the input file, -o is the output file, -q is our minimum quality score and -l is the minimum read length.  The trimmed output files are what we will be using for the next steps of our analysis.
```3. Alignment of Trimmed Reads Using STAR:```

For this next step, you will first need to download the reference genome and annotation file for Glycine max (soybean). The files I used can be found at the following link:

Phytozome – Glycine max
Now that you have the genome and annotation files, you will create a genome index using the following script:

```
STAR --runMode genomeGenerate --genomeDir /common/RNASeq_Workshop/Soybean/gmax_genome/ --genomeFastaFiles /common/RNASeq_Workshop/Soybean/gmax_genome/Gmax_275_v2 --sjdbGTFfile /common/RNASeq_Workshop/Soybean/gmax_genome/Gmax_275_Wm82.a2.v1.gene_exons --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --runThreadN 8
```

The assembly file, annotation file, as well as all of the files created from indexing the genome can be found in ~/datasets

```4. Mapping trimmed reads:```
Now that you have your genome indexed, you can begin mapping your trimmed reads with the following script:The script for mapping all six of our trimmed reads to .bam files can be found in mapped.sh
```5. Convert BAM Files to Raw Counts with HTSeq```: We will use HTSeq to transform these mapped reads into counts that we can analyze with R. “-s” indicates we do not have strand specific counts. “-r” indicates the order that the reads were generated, for us it was by alignment position. “-t” indicates the feature from the annotation file we will be using, which in our case will be exons. The script for generating the count file is htseq_count.sh

```
htseq-count -s no -r pos —t exon -i pacid -f bam SRR391535Aligned.sortedByCoord.out.bam /common/RNASeq_Workshop/Soybean/gmax_genome/Gmax_275_Wm82.a2.v1.gene_exons > SRR391535-output_basename.counts
```

```6. Analysis of Counts with DESeq2```:

For the remaining steps, you can either run the analysis on R in computer or server. You will also need to download R to run DESeq2, and I’d also recommend installing RStudio, which provides a graphical interface that makes working with R scripts much easier.
The R DESeq2 library also must be installed. To install this package, start the R console and enter:

```
source("http://bioconductor.org/biocLite.R")
```

```
biocLite("DESeq2")
```

The R code Deseq_analysis.R do DeSeq2 analysis to get normalize gene count and generate “-replaceoutliers-results.csv” which has adjusted and normal p-values, as well as log2foldchange for all of the genes.

The R code visualization.R  do exploratory data analysis of RNAseq data with DESeq2 do a variety of visualization, QC and other plots to get a sense of what the RNAseq data looks like based on DESEq2 analysis
The plots includes:
 1) MA plot
 2) rlog stabilization and variance stabiliazation
 3) variance stabilization plot
 4) heatmap of clustering analysis
 5) PCA plot

The .csv output file that you get from this R code should look something like this:

![alt text](https://bioinformatics.uconn.edu/wp-content/uploads/sites/15/2015/07/Screen-Shot-2015-07-13-at-1.11.41-PM.png)

Below are some examples of the types of plots you can generate from RNAseq data using DESeq2:
![alt text](https://bioinformatics.media.uconn.edu/wp-content/uploads/sites/15/2015/07/Rplot-logfoldchange.png)
![alt text](https://bioinformatics.uconn.edu/wp-content/uploads/sites/15/2015/07/Gmax_DESeq2-Dispersion-plot.png)
![alt text](https://bioinformatics.uconn.edu/wp-content/uploads/sites/15/2015/07/Rplot-heatmap1.png)
![alt text](https://bioinformatics.uconn.edu/wp-content/uploads/sites/15/2015/07/Rplot-pca.png)

 

7. Merging Data and Using Biomart:

To continue with analysis, we can use the .csv files we generated from the DeSEQ2 analysis and find gene ontology. This next script Biomart.R contains the actual biomaRt calls, and uses the .csv files to search through the Phytozome database. If you are trying to search through other datsets, simply replace the “useMart()”  command with the dataset of your choice. Again, the biomaRt call is relatively simple, and this script is customizable in which values you want to use and retrieve.

After fetching data from the Phytozome database based on the PAC transcript IDs of the genes in our samples, a .txt file is generated that should look something like this:
![alt_text](https://bioinformatics.uconn.edu/wp-content/uploads/sites/15/2015/07/Screen-Shot-2015-07-20-at-6.54.17-PM.png)
