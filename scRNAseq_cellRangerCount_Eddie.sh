
#!/bin/bash
# Grid Engine options (lines prefixed with #$)
# job name: -N
# use the current working directory: -cwd
# number of cores -pe sharedmem
# runtime limit: -l h_rt
# memory limit: -l h_vmem
#$ -N tumourHamedscMouse
#$ -cwd
#$ -pe sharedmem 8
#$ -l h_rt=48:00:00
#$ -l h_vmem=20G
#$ -e error.txt
#$ -o output.txt

## Initialising the environment modules
. /etc/profile.d/modules.sh
# /exports/applications/apps/user-scripts/

# Loading modules*
module load igmm/apps/cellranger/7.1.0/
module load igmm/apps/STAR/2.7.11b

#Running programs
REF_PATH=//exports/eddie/scratch/ngiannak/10x/mouse/refdata-GRCm38-1.0.0

FQ_PATH=//exports/eddie/scratch/ngiannak/Data/HamedNatComms/fastqTumourData
CR_NAME=10032025_HamedTumourscMouse

cellranger count --id=$CR_NAME \
--transcriptome=$REF_PATH \
--fastqs=$FQ_PATH \
--sample=SRR196431,SRR196432 \
