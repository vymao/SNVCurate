#!/bin/bash

#SBATCH -c 1                              # Request core count
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-06:00:00                        # Runtime in D-HH:MM format
#SBATCH -p park                 # Partition to run in
#SBATCH --account=park_contrib
#SBATCH --mem=2000                          # Memory total in MB (for all cores)
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL

module load gcc/6.2.0 samtools

path2SNVCurate="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
bam_dir=$1
output_dir=$2

cd ${output_dir}

for file in ${bam_dir}/*.bam; do
    cd ${output_dir}
	sample=$(basename $file | cut -d'.' -f1)
    if [ ! -f ${output_dir}/${sample}.bam ]; then
          sbatch ${path2SNVCurate}/runRenaming.sh ${file} ${output_dir}
    fi 
done
