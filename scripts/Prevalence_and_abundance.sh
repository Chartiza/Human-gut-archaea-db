#!/bin/bash
#SBATCH -J getAbud_arch
#SBATCH --partition=main
#SBATCH -t 190:00:00
#SBATCH --error=log/coverM_err2
#SBATCH --output=log/coverM.out
#SBATCH --array=0-550%50
#SBATCH --cpus-per-task=24
#SBATCH --mem=120G

module load miniconda3
source activate coverm-singlem-lorikeet

echo 'START: '$(date)

mkdir -p $PWD/abud_results_arch21
output=$PWD/abud_results_arch21

path='Reads/Cleaned_reads_Illumina_2509'
smpls='smpls_lists'
ref='Archaea_ESTrep-21/bins/'

# Find all "reads" files in all subdirectories, excluding the specified subfolders
SUB_LIST=($(<$smpls/illumina_reads_smpls_list_ae))
SUB=${SUB_LIST[${SLURM_ARRAY_TASK_ID}]}
# remove leading and trailing whitespace (spaces or tabs) 
SUB=$(echo "$SUB" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')  
#define sample name
NAME=$(echo "$SUB" | cut -d'/' -f2)

echo ${path}/${SUB}_nohost.fq1.gz
    
coverm genome --coupled ${path}/${SUB}_nohost.fq1.gz ${path}/${SUB}_nohost.fq2.gz --genome-fasta-directory $ref -x fa -o ${output}/abud_${NAME}.tsv --threads 24


echo 'END: '$(date)
echo '---------------'