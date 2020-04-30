#!/bin/bash

#SBATCH -c 1                              # Request core count
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-06:00:00                        # Runtime in D-HH:MM format
#SBATCH -p short                 # Partition to run in
# --account=park_contrib
#SBATCH --mem=1000                          # Memory total in MB (for all cores)
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL


module load gcc/6.2.0 python/3.6.0 java bcftools

path2SNVCurate="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
out=$1
path2Mutect=$2
reference=$3
path2normal=$5
csv=$4

cd $out

rm -f *reads.vcf *somatic.txt *combined.txt


for path in ${out}/*; do 
    [ -d $path ] || continue
    cd ${path}/slurm_submissions
    if ! [ -z "$path2normal" ]; then
        sbatch ${path2SNVCurate}/runAnnotate.sh $path $out $path2Mutect $reference $path2normal $csv $path2SNVCurate
    else
        sbatch ${path2SNVCurate}/runAnnotate.sh $path $out $path2Mutect $reference False $csv $path2SNVCurate
    fi
done



