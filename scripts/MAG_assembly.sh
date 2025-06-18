#!/bin/bash
#SBATCH -J assembly
#SBATCH --partition=main
#SBATCH -t 190:00:00
#SBATCH --error=err2
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G

module load broadwell/megahit/1.2.9
PWD=$(pwd)

while read FQ; do

  NAMEs=`echo "$FQ" | cut -d'/' -f16`
  echo $NAMEs

  megahit -1 $FQ.fq1.gz -2 $FQ.fq2.gz -t 24 -o $PWD/Contigs_DB/assembly_$NAMEs
  rm -r $PWD/Contigs_DB/assembly_$NAMEs/intermediate_contigs

done <smpls.csv

module purge
