#!/bin/bash
#SBATCH -J gtdbtk
#SBATCH --partition=amd
#SBATCH -t 40:00:00
#SBATCH --error=gtdbtk_err
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G

export GTDBTK_DATA_PATH=/gtdb_db/releases/release214

module load py-gtdbtk/2.3.2
module load py-pydantic
echo $GTDBTK_DATA_PATH

PWD=$(pwd)
Bins='Archaea_ESTrep-21/bins/'
mkdir -p $PWD'/gtdb-v214_ESTrep-21'
output=$PWD'/gtdb-v214_ESTrep-21'

echo 'run GTDB-tk'
echo $(date)
echo '...'

gtdbtk classify_wf --genome_dir $Bins --out_dir $output --cpus 64 --scratch_dir $output'/scratch_dir' -x fa --mash_db $output'/mash_db'

echo $(date)
module purge
