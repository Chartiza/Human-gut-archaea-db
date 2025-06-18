#!/bin/bash
#SBATCH -J cm2
#SBATCH --partition=main
#SBATCH -t 190:00:00
#SBATCH --error=cm2_error
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G

module load any/python
source activate checkm2

PWD=$(pwd)

# Define variables
input_dir="/PerSmplsBins_DB_MultiBinning"
file_list="all_bins_list.txt"
batch_prefix="batch_"
batch_size=5000
batch_dir="cm2_all_bins"

# Create file list
find "$input_dir" -name "*.fa" > "$file_list"

# Create batch directory
mkdir -p "$batch_dir"

# Split file list into batches
split -l "$batch_size" "$file_list" "$batch_dir/$batch_prefix"

# Process each batch
for batch_file in "$batch_dir"/"$batch_prefix"*; do
    echo "Processing $batch_file..."
    echo $(date)
    mkdir -p $PWD/$batch_file
    checkm2 predict --threads 34 -x fa --input $(cat "$batch_file") --output-directory $PWD/$batch_file

    echo "---"
done

module purge