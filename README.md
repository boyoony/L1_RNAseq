#!/bin/bash

# 20210519_L1_RNAseq
RNA-sequencing of L1 N2 worms and L1 adr-1(-) worms (Zhao et al., 2015)


##############################
#Downloading SRA files
##############################

#Download SRA files to mac desktop

#Activate rnaseq environment 

conda activate rnaseq

#Unload perl 

module unload perl 

# Variables.
	
ASSEMBLY='ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.genomic.fa.gz'
ANNOTATION='ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/c_elegans/PRJNA13758/c_elegans.PRJNA13758.WS275.canonical_geneset.gtf.gz'


#Make a new directory

cd /N/slate/by6

mkdir L1_RNAseq


#Secure copy and paste files from mac desktop to this folder 
#In the mac terminal, 

scp /Users/yang/Desktop/SRR* by6@carbonate.uits.iu.edu :/N/slate/by6/L1_RNAseq



##############################
##Download SRAtoolkit
##############################

#At your home directory, /N/slate/by6

wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.1/sratoolkit.2.4.1-ubuntu64.tar.gz

tar xzvf sratoolkit.2.4.1-ubuntu64.tar.gz

export PATH=$PATH:/directory/sratoolkit.2.4.1-ubuntu64/bin 


#Convert SRA files to fastq files 

fastq-dump -I --split-files SRR953117.2

fastq-dump -I --split-files SRR953123.2


#Create a directory to store the SRA files

mkdir -p SRAfiles


#Move SRA files to the created folder 

mv SRR953117.2 /N/slate/by6/L1_RNAseq/SRAfiles
mv SRR953123.2 /N/slate/by6/L1_RNAseq/SRAfiles

#Create a directory to store the SRA files

mkdir -p SRAfiles


#Move SRA files to the created folder 

mv SRR953117.2 /N/slate/by6/L1_RNAseq/SRAfiles
mv SRR953123.2 /N/slate/by6/L1_RNAseq/SRAfiles



################################
## Generate STAR Genome Index 
################################
	
# Make a directory to store the genome files in the L1_RNAseq folder 
	
mkdir -p genome

cd /N/slate/by6/L1_RNAseq

	
	# Download and unpack the genome assembly.
	
curl $ASSEMBLY | gunzip > ./genome/assembly.fasta
	
	# Download and unpack the genome annotation.
	
curl $ANNOTATION | gunzip > ./genome/annotation.gtf
	

	# Create a directory to store the index.
	
mkdir -p genome/index

	# Create the STAR genome index.
	# --genomeSAindexNbases 12 was recommended by software.
	
	STAR \
	  --runThreadN 4 \
	  --runMode genomeGenerate \
	  --genomeDir genome/index \
	  --genomeFastaFiles genome/assembly.fasta \
	  --sjdbGTFfile genome/annotation.gtf \
	  --genomeSAindexNbases 12
	


###########################
## Align Reads to Genome ##
###########################
	
	# Create an output directory for aligned reads.
	
mkdir -p results/aligned
	
	# Align the reads.
	
FASTQ=$SRR*
	
	for FASTQ in ${FASTQ[@]}; do
	  PREFIX=results/aligned/$(basename $FASTQ .fastq)_
	  STAR \
	    --runThreadN 8 \
	    --outFilterMultimapNmax 1 \
	    --outFilterScoreMinOverLread .66 \
	    --outFilterMismatchNmax 10 \
	    --outFilterMismatchNoverLmax .3 \
	    --runMode alignReads \
	    --genomeDir genome/index \
	    --readFilesIn $FASTQ \
	    --outFileNamePrefix $PREFIX \
	    --outSAMattributes All \
	    --outSAMtype BAM SortedByCoordinate
	done
	
  
  
############################
#Downloading SAILOR Program
############################


#In Carbonate
#Make directory for analysis

cd /geode2/home/u010/by6/Carbonate

mkdir L1_RNAseq
	

#Load the singularity module

module load singularity


#scp sailor program downloaded from Yeo lab github/HundleyLab project space
#In HundleyLab project space, /N/project/HundleyLab	

scp sailor-1.0.4 by6@carbonate.uits.iu.edu:/geode2/home/u010/by6/Carbonate


############################
#Moving all files needed for SAILOR to Carbonate (not in a directory) 
#############################

#Go to genome directory 

cd /N/slate/by6/L1_RNAseq/genome


#copy and paste assembly.fasta files to Carbonate 

cp assembly.fasta /geode2/home/u010/by6/Carbonate


#Go to L1_RNAseq folder

cd /N/slate/by6/L1_RNAseq


#Make a folder for SAILOR

mkdir sailor 


#Go to aligned directory

cd /N/slate/by6/L1_RNAseq/results/aligned 


# create merged bam file (only if you have multiple replicates) 

module load samtools

samtools merge -@ 8 ../../sailor/SRR953117.2_1_Aligned.sortedByCoord.out.bam SRR953117.2_2_Aligned.sortedByCoord.out.bam

samtools merge -@ 8 ../../sailor/SRR953123.2_1_Aligned.sortedByCoord.out.bam SRR953123.2_2_Aligned.sortedByCoord.out.bam


#Go to Carbonate 

cd /geode2/home/u010/by6/Carbonate


#Copy and paste merged bam files 

cp /N/slate/by6/L1_RNAseq/sailor/* .



#########################
#Get C. elegans reference genome file with SNPs
#########################

#On mac desktop, go to HundleyLab server - protocols - Bioinformatics - Bioinformatics notes - Files for Sailor 

#Copy and paste this to your desktop - c.elegans.WS275.snps.sorted.bed

#In mac desktop terminal, copy the file to Carbonate 

scp /Users/yang/Desktop/c.elegans.WS275.snps.sorted.bed by6@carbonate.uits.iu.edu: 
	

#In Carbonate (not in a directory), you should have bam files, snps.sorted.bed file, assembly.fasta file


############################
#Mark downloaded program as executable within directory
############################

chmod +x sailor-1.0.4
	
./sailor-1.0.4
	
y


#######
#Open yaml file with vim and edit filenames
#######

# open .yaml file
	
vim ce11_example.yaml
	
i
	

# update filenames using arrow keys
# You can put only one file each. If you have multiple files, perform SAILOR multiple times. 
	
# save and quit (ESC - shift+ZZ)
	

# rename .yaml file
	
mv ce11_example.yaml N2.yaml
	

# run sailor
	
./sailor-1.0.4 N2.yaml

