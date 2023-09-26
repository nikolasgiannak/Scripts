featureCounts -a /data10/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o counts_for_DEGs.tsv -T 8 -O mm_Donor_subcutan_fat_Rep1.bam mm_Donor_subcutan_fat_Rep2.bam mm_Donor_subcutan_fat_Rep3.bam mm_TransplantedKO_subcutan_fat_Rep1.bam mm_TransplantedKO_subcutan_fat_Rep2.bam mm_TransplantedKO_subcutan_fat_Rep3.bam

cut -f 1,7,8,9,10,11,12 counts_for_DEGs.tsv > gene_counts_donor_transplKO.txt


featureCounts -a /data10/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o counts_for_DEGs.tsv -T 8 -O mm_Donor_subcutan_fat_Rep1.bam mm_Donor_subcutan_fat_Rep2.bam mm_Donor_subcutan_fat_Rep3.bam mm_Sham_operated_CTR_subcutan_fat_Rep1.bam mm_Sham_operated_CTR_subcutan_fat_Rep2.bam mm_Sham_operated_CTR_subcutan_fat_Rep3.bam 

cut -f 1,7,8,9,10,11,12 counts_for_DEGs.tsv > gene_counts_donor_sham.txt

featureCounts -a /data10/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o counts_for_DEGs.tsv -T 8 -O  mm_Sham_operated_CTR_subcutan_fat_Rep1.bam mm_Sham_operated_CTR_subcutan_fat_Rep2.bam mm_Sham_operated_CTR_subcutan_fat_Rep3.bam mm_TransplantedKO_subcutan_fat_Rep1.bam mm_TransplantedKO_subcutan_fat_Rep2.bam mm_TransplantedKO_subcutan_fat_Rep3.bam

cut -f 1,7,8,9,10,11,12 counts_for_DEGs.tsv > gene_counts_sham_trasnplKO.txtfile,condition,type



file,condition,type
mm_Donor_subcutan_fat_Rep1.bam,donor,single-end
mm_Donor_subcutan_fat_Rep2.bam,donor,single-end
mm_Donor_subcutan_fat_Rep3.bam,donor,single-end
mm_TransplantedKO_subcutan_fat_Rep1.bam,transplKO,single-end
mm_TransplantedKO_subcutan_fat_Rep2.bam,transplKO,single-end
mm_TransplantedKO_subcutan_fat_Rep3.bam,transplKO,single-end


file,condition,type
mm_Donor_subcutan_fat_Rep1.bam,donor,single-end
mm_Donor_subcutan_fat_Rep2.bam,donor,single-end
mm_Donor_subcutan_fat_Rep3.bam,donor,single-end
mm_Sham_operated_CTR_subcutan_fat_Rep1.bam,sham,single-end
mm_Sham_operated_CTR_subcutan_fat_Rep2.bam,sham,single-end
mm_Sham_operated_CTR_subcutan_fat_Rep3.bam,sham,single-end

file,condition,type
mm_Sham_operated_CTR_subcutan_fat_Rep1.bam,sham,single-end
mm_Sham_operated_CTR_subcutan_fat_Rep2.bam,sham,single-end
mm_Sham_operated_CTR_subcutan_fat_Rep3.bam,sham,single-end
mm_TransplantedKO_subcutan_fat_Rep1.bam,transplKO,single-end
mm_TransplantedKO_subcutan_fat_Rep2.bam,transplKO,single-end
mm_TransplantedKO_subcutan_fat_Rep3.bam,transplKO,single-endf


##calculate counts, create need files for DESEQ2 analysis
##HD1_HighD2
featureCounts -a /data10/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o counts_for_DEGs.tsv -T 8 -O mm_Ly6Chigh_D1_1.bam mm_Ly6Chigh_D1_2.bam mm_Ly6Chigh_D1_3.bam mm_Ly6Chigh_D2_1.bam mm_Ly6Chigh_D2_2.bam mm_Ly6Chigh_D2_3.bam
cut -f 1,7,8,9,10,11,12 counts_for_DEGs.tsv > gene_counts_HD1_HighD2.txt

file,condition,type
mm_Ly6Chigh_D1_1.bam,ctrl,single-end
mm_Ly6Chigh_D1_2.bam,ctrl,single-end xmm_Ly6Chigh_D1_3.bam,ctrl,single-end
mm_Ly6Chigh_D2_1.bam,hd2,single-end
mm_Ly6Chigh_D2_2.bam,hd2,single-end
mm_Ly6Chigh_D2_3.bam,hd2,single-end

##HD1_HighD4

featureCounts -a /data10/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o counts_for_DEGs.tsv -T 8 -O mm_Ly6Chigh_D1_1.bam mm_Ly6Chigh_D1_2.bam mm_Ly6Chigh_D1_3.bam mm_Ly6Chigh_D4_1.bam mm_Ly6Chigh_D4_2.bam
cut -f 1,7,8,9,10,11 counts_for_DEGs.tsv > gene_counts_HD1_HighD4.txt

file,condition,type
mm_Ly6Chigh_D1_1.bam,ctrl,single-end
mm_Ly6Chigh_D1_2.bam,ctrl,single-end
mm_Ly6Chigh_D1_3.bam,ctrl,single-end
mm_Ly6Chigh_D4_1.bam,hd4,single-end
mm_Ly6Chigh_D4_2.bam,hd4,single-end

##HD1_LowD2

featureCounts -a /data10/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o counts_for_DEGs.tsv -T 8 -O mm_Ly6Chigh_D1_1.bam mm_Ly6Chigh_D1_2.bam mm_Ly6Chigh_D1_3.bam mm_Ly6Clow_D2_1.bam mm_Ly6Clow_D2_2.bam mm_Ly6Clow_D2_3.bam
cut -f 1,7,8,9,10,11,12 counts_for_DEGs.tsv > gene_counts_HD1_LowD2.txt

file,condition,type
mm_Ly6Chigh_D1_1.bam,ctrl,single-end
mm_Ly6Chigh_D1_2.bam,ctrl,single-end
mm_Ly6Chigh_D1_3.bam,ctrl,single-end
mm_Ly6Clow_D2_1.bam,ld2,single-end
mm_Ly6Clow_D2_2.bam,ld2,single-end
mm_Ly6Clow_D2_3.bam,ld2,single-end

##HD1_LowD4

featureCounts -a /data10/genomes/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf -o counts_for_DEGs.tsv -T 8 -O mm_Ly6Chigh_D1_1.bam mm_Ly6Chigh_D1_2.bam mm_Ly6Chigh_D1_3.bam mm_Ly6Clow_D4_1.bam mm_Ly6Clow_D4_2.bam
cut -f 1,7,8,9,10,11 counts_for_DEGs.tsv > gene_counts_HD1_LowD4.txt

file,condition,type
mm_Ly6Chigh_D1_1.bam,ctrl,single-end
mm_Ly6Chigh_D1_2.bam,ctrl,single-end
mm_Ly6Chigh_D1_3.bam,ctrl,single-end
mm_Ly6Clow_D4_1.bam,ld4,single-end
mm_Ly6Clow_D4_2.bam,ld4,single-end


