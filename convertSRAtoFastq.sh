#Load sra toolkit module and configure it

module load igmm/apps/sratoolkit/2.10.8  

#Run accessions are used to download SRA data. To download a list of Run accessions selected from your Entrez search 
# Click Send to on the top of the page, check the radiobutton File, select Accession List.
#	•	Save this file in the location from which you are running the SRA Toolkit. SraAccList.txt 

#In your command line write

prefetch --option-file SraAccList.txt	

#Convert sra files to fastq files as follows. Add --split-files if you're dealing with paired-end reads.

fasterq-dump  SRR*


OR

prefetch --option-file SraAccList.txt  && fastq-dump SRR*



# VLOOKUP(A3,Neutrophil_data_vsd.txt!$1:$1048576,34,FALSE)