# v1 2021-03-20 based on v 44 2019-05-03

fprimer=GTGYCAGCMGCCGCGGTAA
fprimlength=19
rprimer_rc=
AAACTYAAAKRAATTGRCGG

#Length of amplicon minus primers (908-534 = 374)
crop_length=369

#Make directory where the merged read files will be placed (fastq and fasta)
mkdir merged_reads

#Iterate over all subdirectories, assume they have the sample name, enter and
#uncompress fastq files
for d in $*; do
    echo $d
    cd $d
    pwd
    i=${d//\/}
    echo $i
    
    # Merge the reads (non-staggering since adaptors removed) max 5 diff

    echo $d >> ../readprep.log

    vsearch --fastq_mergepairs *_R1.fastq.gz --reverse *_R2.fastq.gz  --fastq_allowmergestagger --fastqout ${i}_m.fastq --fastq_maxdiffs 20 2>> ../readprep.log

    # 10% errors max (default) and no missing bases
    cutadapt -g $fprimer -m 100 --max-n=0 --discard-untrimmed -O $fprimlength -o ${i}_trimmed1.fastq  ${i}_m.fastq >> ../readprep.log 

    cutadapt -a $rprimer_rc -m 100 --discard-untrimmed -O 19 -o ${i}_trimmed.fastq ${i}_trimmed1.fastq >>../readprep.log

    rm  ${i}_trimmed1.fastq

    # Max expected errors 1
    vsearch --fastq_filter ${i}_trimmed.fastq --fastaout ${i}_merged_QF.fasta --fastq_maxee 1 --fastq_trunclen $crop_length 2>>../readprep.log    
    
    #Fix read names to include "barcode" label
   python ~/script/drive5/fixreads.py ${i}_merged_QF.fasta $i > ../merged_reads/${i}_reads_fixed.fasta

   rm ${i}_merged_QF.fasta

   cd ..
done

echo "Preparing FASTQC reports"

cat */*_m.fastq > all_merged.fastq
fastqc all_merged.fastq
rm all_merged.fastq

cat */*_trimmed.fastq > all_trimmed.fastq
fastqc all_trimmed.fastq
rm all_trimmed.fastq

# cat */*R1*.fastq > all_R1.fastq
# fastqc all_R1.fastq
# rm all_R1.fastq

# cat */*R2*.fastq > all_R2.fastq
# fastqc all_R2.fastq
# rm all_R2.fastq

rm */*.fastq
rm *fastqc.zip

mkdir ../fastQC
mv *.html ../fastQC

for h in ../fastQC/*.html; do gnome-open $h; done

for ff in merged_reads/*.fasta; do
    echo "${ff//_reads_fixed.fasta} `grep -c ">" $ff`" >> readprep.log
done
