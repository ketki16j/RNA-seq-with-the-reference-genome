# RNA-seq-with-the-reference-genome

This repository is a tutorial for differential expression analysis using RNA-Seq data
We will use a subset of gene expression data from Reid et al. 2016, a study of evolutionary adaptation to industrial pollutants in the coastal fish Fundulus heteroclitus. The study consisted of an experimental exposure of embryonic fish from eight populations to the toxicant PCB-126. Here we will use exposed and control samples from two populations.
The samples we will be using are described by the following accession numbers; SRR391535, SRR391536, SRR391537, SRR391538, SRR391539, and SRR391541.

We will be going through quality control of the reads, alignment of the reads to the reference genome, conversion of the files to raw counts, analysis of the counts with DeSeq2, and finally annotation of the reads using Biomart.

1. The script for downloading .SRA files and converting them to fastq can be found in ~/scripts/fastq.sh
2. Quality Control on the Reads Using Sickle: Step one is to perform quality control on the reads using Sickle. We are using unpaired reads, as indicated by the “se” flag in the script quality_control.sh. The -f flag designates the input file, -o is the output file, -q is our minimum quality score and -l is the minimum read length.  The trimmed output files are what we will be using for the next steps of our analysis.
3. Alignment of Trimmed Reads Using STAR:

For this next step, you will first need to download the reference genome and annotation file for Glycine max (soybean). The files I used can be found at the following link:

Phytozome – Glycine max
Now that you have the genome and annotation files, you will create a genome index using the following script:

```
STAR --runMode genomeGenerate --genomeDir /common/RNASeq_Workshop/Soybean/gmax_genome/ --genomeFastaFiles /common/RNASeq_Workshop/Soybean/gmax_genome/Gmax_275_v2 --sjdbGTFfile /common/RNASeq_Workshop/Soybean/gmax_genome/Gmax_275_Wm82.a2.v1.gene_exons --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 100 --runThreadN 8
```

The assembly file, annotation file, as well as all of the files created from indexing the genome can be found in ~/datasets
4. Mapping trimmed reads:
Now that you have your genome indexed, you can begin mapping your trimmed reads with the following script:The script for mapping all six of our trimmed reads to .bam files can be found in mapped.sh
5. Convert BAM Files to Raw Counts with HTSeq: We will use HTSeq to transform these mapped reads into counts that we can analyze with R. “-s” indicates we do not have strand specific counts. “-r” indicates the order that the reads were generated, for us it was by alignment position. “-t” indicates the feature from the annotation file we will be using, which in our case will be exons. The script for generating the count file is htseq_count.sh
