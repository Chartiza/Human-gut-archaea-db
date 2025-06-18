#!/bin/bash
#SBATCH -J BnClv2
#SBATCH --partition=main
#SBATCH -t 40:00:00
#SBATCH --error=BinsCl_errFn
#SBATCH --cpus-per-task=18
#SBATCH --mem=64G

#conda activate default
module load mummer/3.23 any/ANIcalculator/v1

PWD=$(pwd)
mkdir -p clustr
output="${PWD}/clustr"
bins="${PWD}/bins"

echo 'Start clusterins ...'
echo $(date)

dRep dereplicate -p 18 -sa 0.95 --checkM_method taxonomy_wf -comp 25 -con 50 --S_algorithm fastANI --multiround_primary_clustering --greedy_secondary_clustering --run_tertiary_clustering $output -g $bins/*.fa

echo $(date)
module purge
