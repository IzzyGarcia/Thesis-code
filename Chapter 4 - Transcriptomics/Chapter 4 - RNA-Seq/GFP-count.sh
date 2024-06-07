GFP-count using R1 file as an example

#Download GFP fasta file
#HISTA2 Indexing step
hisat2-build -p 20 /home/irg6/GFP_counts/GFP.fasta /home/irg6/GFP_counts/GFP.fasta 

#HISAT2 Aligning step
hisat2 -p 20 -x /home/irg6/GFP_counts/GFP.fasta -1 R1/R1_1.fq.gz -2 R1/R1_2.fq.gz |~/Tools/samtools-1.14/samtools view -b -F 4 -o R1_GFP_mapped.bam

#Sort bam files
for sample in *.bam;do echo -en $sample"\n"; name=$(basename $sample ".bam"); ~/Tools/samtools-1.14/samtools sort -@ 10 $sample -o ${name}_sort.bam; done

#samtools count of GFP 
~/Tools/samtools-1.14/samtools view -c -F 4 R1_GFP_mapped_sort.bam

#Visualise data
~/Tools/samtools-1.14/samtools view -c -F 4 R1_GFP_mapped_sort.bam|less -s
