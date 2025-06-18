#!/bin/bash
#SBATCH -J remove_host
#SBATCH --partition=main
#SBATCH -t 39:00:00
#SBATCH --error=remove_host_err1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

module load fastp/0.23.2

# Create an array with all the filenames in the folder
files=(/reads_without_host/*/*R1.fastq.gz)

for file in "${files[@]}"; do
  bNAME=`echo "$file" | cut -d'/' -f11`
  Rfld=`echo "$file" | cut -d'/' -f1-11`
  R=$(echo "$file" | sed 's/.\{7\}$//')
  R1=$R'1.fq.gz'
  R2=$R'2.fq.gz'
  echo 'Run bowtie2 for '$bNAME
  echo 'start time is '$(date)

  fastp \
  -i $R1 -I i$R1 -o $R1.trimmed_R1.fastq.gz -O $R1.trimmed_R2.fastq.gz 
  --detect_adapter_for_pe \
  --cut_front \
  --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --length_required 50 \
  --thread 4 \
  --html report.html

  echo 'end time is '$(date)

done

module purge