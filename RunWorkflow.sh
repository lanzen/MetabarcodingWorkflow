#!/bin/bash

# Anders Lanzen 2024-06-02 v2
# This script uses vsearch for dereplication of each individual sample
# then merges them using code from SLIM, retaining read mapping, allowing
# for an OTU table to be constructed after SWARM clustering. Singletons are
# removed by default (change minAbundance to 1 to retain) and chimeras removed
# using Silvamod and de novo.  Rmoved OTUs are
# removed from the OTU table.
#
# v2 uses own script (parse_swarm.py) instead of relying on SLIM's
# script (matrix_creation) that does not work with SWARM3 output bc
# the .uc UCLUST output is not ordered in the same way as .fasta
# which makes this script fail. This script also renames ISUs

minAbundance=10 # Change to 1 to retain singletons
db_dir=/share/projects/PhyloRefDB
scriptdir=$HOME/script

#Dereplicate each file separately instead of concatenating at first
for fa in *_fixed.fa*; do 
   sed 's/ 1:N:0:.*;/;/g' $fa > ${fa//.*}_f.fasta
   vsearch -derep_fulllength ${fa//.*}_f.fasta -output ${fa//_reads_fixed.fa*}.dfa -sizeout
   rm ${fa//.*}_f.fasta
done

# # Use fasta_merging script from SLIM to build origins.tsv with all sample
# # abundances for unique sequences (protOTUs)
fasta_merging.py *.dfa
rm *.dfa

# # Remove reads with N (251 bp long Caspian sea seq used as is, and mistake in cutadapt max_n=1
# # for QLD estuary sequences
vsearch --fastx_filter merged.fasta --fastq_maxn 0 --fastaout merged_noN.fasta 

# #Abundance sort and retain singletons 
vsearch -sortbysize merged_noN.fasta -output merged_sorted.fasta -minsize 1
# rm merged.fasta

# # Run SWARM
swarm -f -z -t 16 -w SWARM_OTUs.fasta merged_sorted.fasta -o SWARM.swarms -s SWARM.stats

# # Create OTU table
# # will create SWARM_table.tsv and SWARM_OTUs_f.fasta
# # using files named origins.tsv, SWARM.swarms and SWARM_OTUs.fasta as input
parse_swarm.py
rm merged_sorted.fasta
gzip SWARM_OTUs.fasta

# # Sort by size and remove singleton swarms
 vsearch --sortbysize SWARM_OTUs_f.fasta -output SWARM_OTUs_ns.fasta -minsize $minAbundance --threads 20

# Chimera filtering
vsearch --uchime_ref SWARM_OTUs_ns.fasta -db /share/databases/CREST/ssuome1.1/ssuome1.1.fasta -strand plus -nonchimeras SWARM_OTUs_uchime_ref.fa --threads 20

vsearch --uchime_denovo SWARM_OTUs_uchime_ref.fa --nonchimeras SWARM_OTUs_uchime_ref_denovo.fa --threads 20

awk 'NR==1 {print ; next} {printf /^>/ ? "\n"$0"\n" : $1} END {printf "\n"}' SWARM_OTUs_uchime_ref_denovo.fa  > SWARM_OTUs_L.fasta

sed 's/;size=.*//' SWARM_OTUs_L.fasta > SWARM_OTUs_clean.fasta
#rm SWARM_OTUs_L.fasta

filter_OTUtable_by_seq.py SWARM_table.tsv SWARM_OTUs_clean.fasta > SWARM_table_clean.tsv
gzip SWARM_table.tsv

## Skip LULU (uncomment to include)

#makeblastdb -in SWARM_OTUs_clean.fasta -dbtype nucl
#blastn -db SWARM_OTUs_clean.fasta -outfmt '6 qseqid sseqid pident' -out LULU_match_list.txt -qcov_hsp_perc 95 -perc_identity 95 -query SWARM_OTUs_clean.fasta -num_threads 12
#Rscript $scriptdir/LULU.R
# Classify LULU filtered OTUs using CREST
filter_seq_by_OTUtable.py SWARM_OTUs_clean.fasta SWARM_table_curated.tsv > SWARM_OTUs_curated.fasta

# Classify with CREST (v4 and locally or moleco v3 - comment)

conda init bash
conda activate crest
crest4 --fasta SWARM_OTUs_clean.fasta --search_db ssuome -t 20 --otu_table SWARM_table_clean.tsv 

# delete some old files
#tar czf readlinks.tar.gz *_reads_fixed.fasta
#rm *_reads_fixed.fasta





