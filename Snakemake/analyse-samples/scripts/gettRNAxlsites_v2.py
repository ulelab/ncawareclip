# Charlotte Capitanchik
# 25.04.18
# Based on fastCLIP tRNA strategy
# input is SAM of tRNA mapping
# and .fai of what is was mapped to


import pysam
from collections import defaultdict
import operator
from collections import Counter
import sys
import statistics
import random
import re


# input #
#trna_fai = sys.argv[1] # "/Volumes/lab-luscomben/working/charlotte/miCLIP_UN1710_CS2_delta2/summarising_trna_mapping_style/permanent_files/UCSC_hg38_tRNAs_mature_PTM-deduplicated.fa.fai" # #"/Users/capitac/PhD/trna_index/hg38-mature-tRNAs-deduplicated.fa.fai"
reads = sys.argv[1] #"../../results/maturetRNA_mapped/DDX3_WT_iCLIP_rep2.Aligned.out.sorted.deduplicated.bam"  # #/Volumes/lab-luscomben/working/charlotte/miCLIP_UN1710_CS2_delta2/summarising_trna_mapping_style/results/m1a_786o_noq_15min_20180413.Aligned.out.sorted.deduplicated.bam" # #"/Users/capitac/PhD/trna_index/m1a_HEK293_beta_1.Aligned.out.sorted.bam"
dir_and_samp = sys.argv[2] #"testing_DDX3_WT_iCLIP_rep2" #  #"/Volumes/lab-luscomben/working/charlotte/testing" #  #"/Users/capitac/PhD/trna_index/m1a_HEK293_beta_1"
samp = sys.argv[3] #"testing_DDX3_WT_iCLIP_rep2" #
seed = sys.argv[4] #1 #
fraction_merge = float(sys.argv[5]) #0.9 #

# output #
# "_trna_summary_stats.tsv"
# "_trna_mixedup_anticodons.tsv"
# "_ambig_ac_position_dist.tsv"

## Code to parse the .fai file to get the name of all the amino acids and counts .etc ##
## Not required for main output
# trna_index = open(trna_fai).readlines()
# aminoacid_count = defaultdict(int)
# anticodon_list=[]
# aminoacid_list=[]
# trna_index
# for i in range(0,len(trna_index)):
#     linelist = trna_index[i].split("\t")
#     namelist = trna_index[i].split("-")
#     anticodon = "-".join(namelist[1:3])
#     amino_acid = namelist[1]
#     anticodon_list.append(anticodon)
#     aminoacid_list.append(amino_acid)
#     aminoacid_count[amino_acid] += 1
# anticodon_list = set(anticodon_list)
# aminoacid_list = set(aminoacid_list)
# # Find out the amino acid with the highest number of genes
# max_key = max(aminoacid_count.iterkeys(), key=(lambda key: aminoacid_count[key]))
# aminoacid_count[max_key]
#print(aminoacid_count[max_key])

## Now move on to the reads ##

samfile = pysam.AlignmentFile(reads,"rb")

## Summarise at three levels
read_dictionary_amino_acid = defaultdict(set)
read_dictionary_anticodon = defaultdict(set)
read_dictionary_start_positions = defaultdict(list)
read_gene_dict = defaultdict(list)
readname_list = []

for read in samfile.fetch():
    if read.is_unmapped == False:
        aa_tmp = samfile.get_reference_name(read.reference_id)
        aa_tmp = re.sub("^(nm-|nmt-)","", aa_tmp).split("-")
        aa = aa_tmp[1]
        anticodon = "-".join(aa_tmp[1:3])
        read_dictionary_amino_acid[read.query_name].add(aa)
        read_dictionary_anticodon[read.query_name].add(anticodon)
        read_dictionary_start_positions[read.query_name] += [read.reference_start]
        read_gene_dict[read.query_name] += [samfile.get_reference_name(read.reference_id)]
        readname_list += [read.query_name]
    else:
        continue

## #### For a file where anticodons are not merged ######
FINAL_read_dictionary_notmerge = defaultdict(set)
for key, value in read_dictionary_anticodon.items(): 
    if len(value) == 1:
        FINAL_read_dictionary_notmerge[key] = value

final_cdna_counts_dict_notmerge = defaultdict(list)
random.seed(seed)
# assign positions, use mode if possible, but if no mode then select randomly
# from multimapped positions using random seed
for key,value in FINAL_read_dictionary_notmerge.items():
    try:
        x = statistics.mode(read_dictionary_start_positions[key])
        final_cdna_counts_dict_notmerge["".join(value)] += [x]
    except:
        final_cdna_counts_dict_notmerge["".join(value)] += [random.choice(read_dictionary_start_positions[key])]

for key,value in final_cdna_counts_dict_notmerge.items():
    final_cdna_counts_dict_notmerge[key] = Counter(value)

########### Unmerged ends here #################################
## Find out which anticodons are ambiguous ## 
number_of_reads_multimapped_anticodon_dict = defaultdict(int)

for key, value in read_dictionary_anticodon.items():
    if len(value) > 1:
        number_of_reads_multimapped_anticodon_dict[";".join(sorted(value))] += 1

sorted_mixup_ranking = sorted(number_of_reads_multimapped_anticodon_dict.items(), key=operator.itemgetter(1), reverse=True)
total = sum(number_of_reads_multimapped_anticodon_dict.values())
fract_total = total * fraction_merge
print(total)
print(fract_total)

counting=0
thresholded_multimapped_anticodon_dict = defaultdict(int)

# Threshold mixed up anticodons based on % reads we want to recover
for i in sorted_mixup_ranking:
    counting += i[1]
    print(i)
    name=i[0]
    number=i[1]
    thresholded_multimapped_anticodon_dict[name] = number
    if counting >= fract_total:
        break

# This is the actual fraction of reads we are recovering
real_frac = float(sum(thresholded_multimapped_anticodon_dict.values())) / float(total)

# The anticodons we are going to merge
comparison_list = list(thresholded_multimapped_anticodon_dict.keys())

## Merge anticodons that are ambiguous################
FINAL_read_dictionary = defaultdict(set)

# If "readmap" has components in "checklist", then replace with the
# corresponsing entry in checklist ->
# If it doesn't - return readmap
# If it has conflicting entries - return readmap

def getIndex(readmap, checklist):
    """
    inputs:
        readmap - list of symbols
        checklist - list of strings(delineated by semicolons)
    """
    readmap = list(readmap) #type cast from set to list
    holder = []
    counter=0
    
    for i in readmap: # iterate over readmap
        switch=True
        for idx, j in enumerate(checklist):
            k = j.split(';')
            if i in k:
                holder.append(idx)
                switch=False
        if switch is True:
            holder.append(-1)

    if len(holder) < 1:
        holder = [-1]
    #evaluate holder values and return based on this
    holder = set(holder)
    holder = list(holder)

    if -1 in holder:
        return readmap
    elif len(holder) > 1:
        return readmap
    elif len(holder) == 1:
        return [checklist[holder[0]]] #checklist[list(holder)[0]]
    else:
        print('something went wrong')


tempdict = {}

for key, value in read_dictionary_anticodon.items():
    newval = getIndex(value, comparison_list)
    tempdict[key] = newval

# Filter for single mapping anticodon/merged anticodons
for key, value in tempdict.items(): 
    if len(value) == 1:
        FINAL_read_dictionary[key] = value

# assign possitions
final_cdna_counts_dict = defaultdict(list)
random.seed(seed)
for key,value in FINAL_read_dictionary.items():
    try:
        x = statistics.mode(read_dictionary_start_positions[key])
        final_cdna_counts_dict["".join(value)] += [x]
    except:
        final_cdna_counts_dict["".join(value)] += [random.choice(read_dictionary_start_positions[key])]

# count at positions
for key,value in final_cdna_counts_dict.items():
    final_cdna_counts_dict[key] = Counter(value)

# Write a file condensed to groups of anticodons
finalfile = open(dir_and_samp + "_tRNA.bed", "w")
for key,value in final_cdna_counts_dict.items():
    gene = key
    for k in value:
        finalfile.write(gene+"\t"+str(k)+"\t"+str(k+1)+"\t"+gene+"\t"+str(value[k])+"\t+\n")
finalfile.close()

# Write a file where only reads unambiguous for anticodon are used

finalfile = open(dir_and_samp + "_tRNA_unambig_AC.bed", "w")
for key,value in final_cdna_counts_dict_notmerge.items():
    gene = key
    for k in value:
        finalfile.write(gene+"\t"+str(k)+"\t"+str(k+1)+"\t"+gene+"\t"+str(value[k])+"\t+\n")
finalfile.close()

##### ESSENTIAL CODE STOPS HERE: JUST SUMMARY FILES FROM NOW ON  #####

# Some statistics for a summary file

total_mapped_reads = len(read_dictionary_amino_acid)
total_uniquely_mapped = 0
number_of_reads_multimapped_aa = 0 # number of reads that are ambigous for amino acid
number_of_reads_multimapped_anticodon = 0 # number of reads that are ambiguous for anticodon
number_of_reads_multimapped_anticodon_dict = defaultdict(int)

for key, value in read_dictionary_start_positions.items():
    if len(value) == 1:
        total_uniquely_mapped += 1

for key, value in read_dictionary_amino_acid.items():
    if len(value) > 1:
        number_of_reads_multimapped_aa += 1

reads_unambig_AC = defaultdict(str)
reads_unambig_pos = defaultdict(list)
for key, value in read_dictionary_anticodon.items():
    if len(value) > 1:
        number_of_reads_multimapped_anticodon += 1
        number_of_reads_multimapped_anticodon_dict[";".join(value)] += 1
    elif len(value) == 1:
        reads_unambig_AC[key] = "".join(value)
        reads_unambig_pos[key] = read_dictionary_start_positions[key]


# Sort the dictionary by highest to lowest
mixed_up_anticodons = open(dir_and_samp+"_tRNA_mixedup_anticodons.tsv", "w")
sorted_mixup_ranking = sorted(number_of_reads_multimapped_anticodon_dict.items(), key=operator.itemgetter(1), reverse=True)
for i in sorted_mixup_ranking:
    if i[1] > 200:
        mixed_up_anticodons.write(str(i[0])+"\t"+str(i[1])+"\t"+str(total_mapped_reads)+"\t"+samp+"\n")
mixed_up_anticodons.close()

summary_stats = open(dir_and_samp+"_tRNA_summary_stats.tsv", "w")
summary_stats.write("Total mapped reads to trna\t"+str(total_mapped_reads)+"\t"+samp+"\n")
summary_stats.write("Total single mapping reads to trna\t"+str(total_uniquely_mapped)+"\t"+samp+"\n")
summary_stats.write("Number of reads that are ambiguous at amino acid level\t"+str(number_of_reads_multimapped_aa)+"\t"+samp+"\n")
summary_stats.write("Percentage of reads used for amino acid summary\t"+str(round(float(total_mapped_reads-number_of_reads_multimapped_aa)/total_mapped_reads,3))+"\t"+samp+"\n")
summary_stats.write("Number of reads that are ambiguous at anticodon level\t"+str(number_of_reads_multimapped_anticodon)+"\t"+samp+"\n")
summary_stats.write("Percentage of reads used for anticodon summary\t"+str(round(float(total_mapped_reads-number_of_reads_multimapped_anticodon)/total_mapped_reads,3))+"\t"+samp+"\n")
summary_stats.write("The actual fraction of ambiguous reads that are recovered by merging anticodons\t"+str(real_frac)+"\n")
summary_stats.write("The anticodons that were merged are\t"+("|").join(comparison_list))

### AAlevel summary
aa_counts = defaultdict(int)
for key, value in read_dictionary_amino_acid.items():
    if len(value) == 1:
        aa_counts["".join(value)] += 1
### Write to file
aa_counts

### Anticodon level summary & also get names of reads that are unique to one anticodon
anticodon_counts = defaultdict(int)
unique_anticodon_read_names = set()
for key, value in read_dictionary_anticodon.items():
    if len(value) == 1:
        anticodon_counts["".join(value)] += 1
        unique_anticodon_read_names.add(key)
anticodon_counts


### How many reads that multimap within one anticodon have different starting positions? What is the maximum distance between them?

# list of reads that have ambiguous starting positions
# (here we assume that if there is a mode the starting position is "unambiguous")
# ie. if positions are 58 58 58 58 58 59 then we would take 58 as the "unambiguous" start

ambig_pos = 0
ambig_pos_reads = []
ambig_distances = []
for key,value in reads_unambig_pos.items():
    if len(set(value)) != 1:
        try:
            x = statistics.mode(value)
        except:
            ambig_pos += 1
            ambig_pos_reads.append(key)
            ambig_distances.append(max(value)-min(value))

anticodons_ambig_pos = dict((k, read_dictionary_anticodon[k]) for k in ambig_pos_reads if k in read_dictionary_anticodon)
ambig_pos_list = []
for key,value in anticodons_ambig_pos.items():
    ambig_pos_list.append("".join(value))

ac_affected = Counter(ambig_pos_list)

summary_stats.write("Number of reads that have an ambiguous position at anticodon level\t"+str(ambig_pos)+"\t"+samp+"\n")
if len(ambig_distances) != 0:
    summary_stats.write("Median distance between ambiguous positions at anticodon level\t"+str(statistics.median(ambig_distances))+"\t"+samp+"\n")
else:
    summary_stats.write("Median distance between ambiguous positions at anticodon level\t"+"NA"+"\t"+samp+"\n")
if len(ambig_distances) != 0:
    summary_stats.write("Biggest distance between ambiguous positions at anticodon level\t"+str(max(ambig_distances))+"\t"+samp+"\n")
else:
    summary_stats.write("Biggest distance between ambiguous positions at anticodon level\t"+"NA"+"\t"+samp+"\n")
summary_stats.close()

ambig_stats = open(dir_and_samp+"_tRNA_ambig_ac_position_dist.tsv", "w")
for key,value in ac_affected.items():
    ambig_stats.write(key+"\t"+str(value)+"\t"+str(ambig_pos)+"\n")
ambig_stats.close()