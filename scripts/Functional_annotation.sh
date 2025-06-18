#!/bin/bash
#SBATCH -J prokka
#SBATCH --partition=main
#SBATCH -t 36:00:00
#SBATCH --error=prokka_err
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G

module load prokka/1.14.6
module load java/11.0.2

PWD=$(pwd)
Bins=$PWD/bins
mkdir -p $PWD/results_prokka

while read FQ; do

  NAME="${FQ:0:-3}"
  echo $NAME
  prokka --dbdir /gpfs/space/home/pantiukh/your/db/loc --kingdom Archaea --prefix $NAME --outdir $PWD/results_prokka/$NAME $Bins/$FQ
  

done </$PWD/archaea_bins.csv
