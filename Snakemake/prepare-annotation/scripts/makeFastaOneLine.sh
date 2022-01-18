#!/bin/bash

# Make all fastas in a directory into one-line versions RECURSIVELY

files=$(ls -d * */* | grep fa)
file_stem_array=()

for file in ${files}; do
    file_stem=$(basename ${file} .fa)
    file_dir=$(dirname ${file})
    file_stem_array+=( "${file_dir}/${file_stem}" )
done

for file in "${file_stem_array[@]}"; do
    echo ${file}
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${file}.fa | tail -n +2 > ${file}.oneline.fa
done

echo Completed
