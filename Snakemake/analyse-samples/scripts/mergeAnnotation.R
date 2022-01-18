#!/usr/bin/env Rscript

# Charlotte Capitanchik 2018 #
# Merge annotation for miCLIP #

library(data.table)
library(dplyr)
library(GenomicRanges)
options(scipen = 999)
options("scipen"=100, "digits"=4)
args = commandArgs(trailingOnly=TRUE) # get command line arguments in a list [args]

# Run command: Rscript --vanilla scripts/mergeAnnotation.R {input.genome} {input.all_counts} {input.small_rna_annot} {input.rDNAannot} {output.finalannotation}

# INPUT #
in_genome_file <- args[1] #"../results/annotation_xlsites/m1a_HEK293_beta_1.genome_overlaps" #
in_final_file <- args[2] #"../results/annotation_xlsites/m1a_HEK293_beta_1.finalcounts"#args[3]
small_rna_lengths <- args[3] #"../permanent_files/small_rna_annotation_final.txt" #args[4]
rDNA_annotation <- args[4] #"../permanent_files/rDNA_annotation_final.txt" #args[5]
output_file <- args[5] #"../results/annotation_xlsites/m1a_HEK293_beta_1.TEST.final_annotation" #args[6]


genome <- fread(in_genome_file)
all <- fread(in_final_file)

nrow(genome) == nrow(all %>% filter(!grepl("snRNA",V1) & !grepl("rRNA",V1) & !grepl("rDNA",V1) & !grepl("-",V1)) %>% select(V1))
small_RNA <- all %>% filter(grepl("snRNA",V1) | grepl("rRNA",V1) | grepl("rDNA",V1) | grepl("-",V1)) 

# Ensure everything is the same format

genome_polished <- data.frame(genome$V1, genome$V2, genome$V3, genome$V21, genome$V5, genome$V6, genome$V15, genome$V21, genome$V21, genome$V10, genome$V16-1, genome$V17, genome$V8,genome$V9)
colnames(genome_polished) <- c("chr", "start", "stop","name", "count","strand","region","gene_id","subtype","repeat_overlap","gene_start", "gene_stop", "repeat_start", "repeat_stop")
genome_polished$name <- gsub('gene_id ".*";gene_name "',"",genome_polished$name)
genome_polished$name <- gsub('";biotype ".*";',"",genome_polished$name)
genome_polished$name <- gsub('None',".",genome_polished$name)

genome_polished$gene_id <- gsub('gene_id "',"",genome_polished$gene_id)
genome_polished$gene_id <- gsub('";gene_name ".*',"",genome_polished$gene_id)

genome_polished$subtype <- gsub('biotype ""','biotype "."',genome_polished$subtype)
genome_polished$subtype <- gsub('gene_id ".*";gene_name ".*";biotype "',"",genome_polished$subtype)
genome_polished$subtype <- gsub('";',"",genome_polished$subtype)

# Separate out the rDNA to annotate each gene in the cluster properly, then merge back #
rDNA <- small_RNA %>% filter(grepl("rDNA", V1))
rest <- small_RNA %>% filter(!grepl("rDNA", V1))

rDNA_a <- fread(rDNA_annotation)

# change to 1-based coordinates for GRanges
rDNA$V2 <- rDNA$V2+1
colnames(rDNA_a) <- c("chr","start","stop","name","score","strand")
colnames(rDNA) <- c("chr","start","stop","name","score","strand")

rDNA_a <- makeGRangesFromDataFrame(rDNA_a,keep.extra.columns=TRUE)
rDNA <- makeGRangesFromDataFrame(rDNA,keep.extra.columns=TRUE)

hits <- findOverlaps(rDNA, rDNA_a)
rDNA$name <- as.character(NA)
rDNA[queryHits(hits)]$name <- rDNA_a[subjectHits(hits)]$name

rDNA <- data.frame(V1=seqnames(rDNA),
                 V2=start(rDNA)-1,
                 V3=end(rDNA),
                 V4=rDNA$name,
                 V5=rDNA$score,
                 V6=strand(rDNA))

rest$V4 <- rest$V1
small_RNA <- rbind(rDNA, rest)

small_RNA_polished <- data.frame(small_RNA$V1, small_RNA$V2, small_RNA$V3, small_RNA$V4, small_RNA$V5, small_RNA$V6, "ncRNA",".","cytoplasmic tRNA",".",".",".",".",".", stringsAsFactors = FALSE)
colnames(small_RNA_polished) <- c("chr", "start", "stop","name", "count","strand","region","gene_id","subtype","repeat_overlap","gene_start", "gene_stop", "repeat_start", "repeat_stop")
small_RNA_polished$subtype[grepl("snRNA",small_RNA_polished$name)] <- "snRNA"
small_RNA_polished$subtype[grepl("rRNA",small_RNA_polished$chr)] <- "rRNA"
small_RNA_polished$subtype[grepl("rDNA",small_RNA_polished$chr)] <- "rRNA"

lengths <- fread(small_rna_lengths,stringsAsFactors = FALSE)
tmp <- left_join(small_RNA_polished, lengths, by=c("name"="V1"))
tmp$gene_start <- tmp$V2
tmp$gene_stop <- tmp$V3
small_RNA_polished <- tmp %>% select(-V2, -V3)

# Mega-merge
final_annot <- rbind(genome_polished, small_RNA_polished)
final_annot[final_annot==-1] <- "."

final_annot$subtype <- gsub(".","intergenic",final_annot$subtype,fixed=TRUE)
final_annot[grepl("IMMATURE",final_annot$chr),]$subtype <- "immature cytoplasmic tRNA"

fwrite(final_annot, file=output_file, quote = FALSE, col.names = FALSE, sep="\t")