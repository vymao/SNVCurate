#!/bin/bash

#SBATCH -c 1                              # Request core count
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-06:00:00                        # Runtime in D-HH:MM format
#SBATCH -p short                 # Partition to run in
# --account=park_contrib
#SBATCH --mem=1000                          # Memory total in MB (for all cores)
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=victor_mao@hms.harvard.edu   # Email to which notifications will be sent

module load gcc/6.2.0 samtools

bam_dir=$1
output_dir=$2

for file in ${bam_dir}/*.bam; do
	sample=$(basename $file | cut -d'.' -f1)
	samtools view -H $file > ${bam_dir}/${sample}_header.sam
	sed -e "s/\tSM:[^[:space:]]*\t/\tSM:${sample}\t/" ${bam_dir}/${sample}_header.sam
	samtools reheader ${bam_dir}/${sample}_header.sam ${file} > ${output_dir}/${sample}.bam
done