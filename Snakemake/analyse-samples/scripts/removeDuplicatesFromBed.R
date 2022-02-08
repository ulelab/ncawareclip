#!/usr/bin/env Rscript

# Charlotte Capitanchik 2022 #
# Remove dupes from bed file #

library(data.table)
library(dplyr)
options(scipen = 999)
options("scipen"=100, "digits"=4)
args = commandArgs(trailingOnly=TRUE) # get command line arguments in a list [args]

# Run command: Rscript --vanilla scripts/removeDuplicatesFromBed.R {input.bed} {output.bed}

# INPUT #
input_file <- args[1]
output_file <- args[2]


input_bed <- fread(input_file) 

out_fwd = input_bed %>%
                filter(V6=="+") %>%
                mutate(rbc=gsub(".*rbc:","",V4))
                dplyr::group_by(V1, V2, rbc) %>%
                dplyr::slice(n=1)

out_rev = input_bed %>%
                filter(V6=="-") %>%
                mutate(rbc=gsub(".*rbc:","",V4))
                dplyr::group_by(V1, V3, rbc) %>%
                dplyr::slice(n=1)

out = rbind(out_fwd, out_rev)

write.table(out, file=output_file, col.names=FALSE, sep="\t", row.names=FALSE, quote=FALSE)
