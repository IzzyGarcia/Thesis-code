#Create an index file
samtools index input.bam

#count per chromosome
samtools idxstats input.bam| awk '{print $1" "$3}' > Chromosome_count_samplename.txt

#or loop
for sample in *.bam;do echo -en $sample"\n"; name=$(basename $sample ".bam"); ~/Tools/samtools-1.14/samtools index -@ 10 $sample; ~/Tools/samtools-1.14/samtools idxstats $sample | awk '{print $1" "$3}' >Chromosome_stats_${name}.txt; done
