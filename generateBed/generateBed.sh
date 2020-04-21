#!/bin/bash

#SBATCH -c 1                              # Request core count
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-06:00:00                        # Runtime in D-HH:MM format
#SBATCH -p short                 # Partition to run in
# --account=park_contrib
#SBATCH --mem=1000                          # Memory total in MB (for all cores)
#SBATCH --mail-type=FAIL           


module load gcc bedtools bedops

path2SNVCurate=$(dirname "$(readlink -f "$0")")
input_bam=$1
reference=$2
csv=$3
out=$4
cut=$5

dir=$(dirname $input_bam)
sample=$(basename $input_bam | cut -d'.' -f1)

bedtools genomecov -bg -ibam ${input_bam} -g ${reference} > ${out}/${sample}.coverage

python3 ${path2SNVCurate}/analyzeCoverage.py -input_dir ${out} -csv ${csv} -cut ${cut} -out ${out}

bedops --everything ${out}/*cov_cut_mergedadjacent.bed \
    | bedmap --count --echo --delim '\t' - \
    | uniq \
    | awk -v OVERLAPS=2 '$1 >= OVERLAPS' \
    | cut -f2- \
    > ${out}/common.bed

bedops --partition ${out}/common.bed | bedmap --count --echo --delim '\t' - ${out}/common.bed | awk '$1 >= 2' | cut -f2- > ${out}/overlaps.bed