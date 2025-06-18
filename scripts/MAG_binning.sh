#!/bin/bash
#SBATCH -J getBins
#SBATCH --partition=main
#SBATCH -t 190:00:00
#SBATCH --error=binning_err1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G

module load minimap2/2.14
module load samtools/1.3.1
module load metabat/2.15
module load maxbin/2.2.7
module load py-vamb/3.0.7
module load das-tool/1.1.4

# Read values from the command-line argument
subdir=$1
NAMEs=$4
R1=$2
R2=$3

Cont='Contigs_DB/assembly_'
runFld=$(pwd)
mkdir -p "$runFld/Bins_DB_MultiBinning"
PWD="$runFld/Bins_DB_MultiBinning"

echo '--------------------'
echo $NAMEs
echo $(date)
echo '...'
echo $Cont$NAMEs/final.contigs.fa 
echo $R1
echo $R2

# Map reads to contigs
echo '...'
echo 'Map reads to contigs'
mkdir -p $PWD/$NAMEs
mkdir -p $PWD/$NAMEs/coverage

minimap2 -d $PWD/$NAMEs/coverage/catalogue.mmi $Cont$NAMEs/final.contigs.fa
minimap2 -t 12 -N 5 -ax sr $PWD/$NAMEs/coverage/catalogue.mmi $R1 $R2 | samtools view -F 4 -b --threads 12 > $PWD/$NAMEs/coverage/$NAMEs.bam
samtools index $PWD/$NAMEs/coverage/$NAMEs.bam $PWD/$NAMEs/coverage/$NAMEs.bam.bai
samtools sort -o $PWD/$NAMEs/coverage/$NAMEs.sorted.bam $PWD/$NAMEs/coverage/$NAMEs.bam
samtools index $PWD/$NAMEs/coverage/$NAMEs.sorted.bam $PWD/$NAMEs/coverage/$NAMEs.sorted.bam.bai

# Get Bins with Metabat2
echo '...'
echo 'run Metabat2'
mkdir -p $PWD/$NAMEs/metabat
jgi_summarize_bam_contig_depths --outputDepth $PWD/$NAMEs/coverage/$NAMEs.metabat.depth.txt $PWD/$NAMEs/coverage/$NAMEs.sorted.bam
metabat2 -i $Cont$NAMEs/final.contigs.fa -a $PWD/$NAMEs/coverage/$NAMEs.metabat.depth.txt -o $PWD/$NAMEs/metabat/metabat -m 1500

# Get Bins with Vamb
echo '...'
echo 'run Vamb'
vamb --outdir $PWD/$NAMEs/vamb/ --fasta $Cont$NAMEs/final.contigs.fa --bamfiles $PWD/$NAMEs/coverage/$NAMEs.bam --minfasta 50000
for f in $PWD/$NAMEs/vamb/bins/* ; do
  bins=$(echo "$f" | rev | cut -d'/' -f1 | rev)
  mv -- "$PWD/$NAMEs/vamb/bins/$bins" "$PWD/$NAMEs/vamb/bins/vamb.$bins"
done

# Get Bins with Maxbin2
echo '...'
echo 'run Maxbin2'
mkdir -p $PWD/$NAMEs/maxbin
python $runFld/depth_for_maxbin.py $NAMEs $runFld
run_MaxBin.pl -thread 24 -min_contig_length 500 -contig $Cont$NAMEs/final.contigs.fa -out $PWD/$NAMEs/maxbin/maxbin -abund $PWD/$NAMEs/coverage/$NAMEs.maxbin.depth.txt

# Choose best bins with a DAS-tool
echo '...'
echo 'run DAS-tool'
mkdir -p $PWD/$NAMEs/dastool
sh $runFld/Fasta_to_Contig2Bin.sh -i $PWD/$NAMEs/metabat/ -e fa > $PWD/$NAMEs/dastool/metabat.scaffolds2bin.tsv
sh $runFld/Fasta_to_Contig2Bin.sh -i $PWD/$NAMEs/vamb/bins/ -e fna > $PWD/$NAMEs/dastool/vamb.scaffolds2bin.tsv
sh $runFld/Fasta_to_Contig2Bin.sh -i $PWD/$NAMEs/maxbin/ -e fasta > $PWD/$NAMEs/dastool/maxbin.scaffolds2bin.tsv
awk 'BEGIN {OFS = "\t"}{print $1,$5}' $PWD/$NAMEs/dastool/vamb.scaffolds2bin.tsv > $PWD/$NAMEs/dastool/vamb1.scaffolds2bin.tsv

DAS_Tool -i $PWD/$NAMEs/dastool/metabat.scaffolds2bin.tsv,$PWD/$NAMEs/dastool/vamb1.scaffolds2bin.tsv,$PWD/$NAMEs/dastool/maxbin.scaffolds2bin.tsv -l metabat,vamb,maxbin -c $Cont$NAMEs/final.contigs.fa -o $PWD/$NAMEs/dastool/bins --write_bins --threads 12


echo $(date)
module purge