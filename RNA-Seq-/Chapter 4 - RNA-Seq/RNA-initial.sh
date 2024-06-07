RNA-Seq analysis using R1 sample

#Downloading data files from Novogene
mkdir RNASeq-data
cd RNASeq-data
wget -c https://englandaws-data1.s3-eu-west-1.amazonaws.com/out/CP2019022700023/X204SC21110661-Z01-F001/MD5.txt
wget -c https://englandaws-data1.s3-eu-west-1.amazonaws.com/out/CP2019022700023/X204SC21110661-Z01-F001/checkSize.xls
wget -c https://englandaws-data1.s3-eu-west-1.amazonaws.com/out/CP2019022700023/X204SC21110661-Z01-F001/X204SC21110661-Z01-F001.zip
md5sum -c MD5.txt >Check_output.txt
unzip X204SC21110661-Z01-F001.zip

#Downloading packages 
python /opt/icarus/setup.py
conda install -c bioconda fastqc 
fastqc -h 

#FastQC of RNA data files 
fastqc R1/R1_1.fq.gz

#Perparing the reference genome (release 105) & downloading Ensemble toplevel fasta file and chromosome GTF file 
conda install -c bioconda hisat2
conda install -c bioconda samtools
wget -c http://ftp.ensembl.org/pub/release-105/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
gunzip Homo_sapiens.GRCh38.dna.toplevel.fa.gz
wget -c http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/Homo_sapiens.GRCh38.105.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.105.chr.gtf.gz

#Using HISAT2 to generate index
hisat2-build -p 10 Homo_sapiens.GRCh38.dna.toplevel.fa Homo_sapiens.GRCh38.dna.toplevel

#Using HISAT2 to map to genome
hisat2 -x /home/irg6/Homo_sapiens.GRCh38.dna.toplevel -1 R1/R1_1.fq.gz -2 R1/R1_2.fq.gz |~/Tools/samtools-1.14/samtools view -b -o R1.bam

#Sort the bam files generated from mapping step (looped)
for sample in *.bam;do echo -en $sample"\n"; name=$(basename $sample ".bam"); ~/Tools/samtools-1.14/samtools sort -@ 10 $sample -o ${name}_sort.bam; done

#Feature counts (pair-end reads)
conda install -c bioconda subread 
featureCounts -p -T5 -t exon -g gene_id -a /home/irg6/Homo_sapiens.GRCh38.105.chr.gtf -o counts_PE.txt R1_sort.bam 

Export counts matric to PC and open in excel format 
