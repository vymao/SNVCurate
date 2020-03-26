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


module load gcc/6.2.0 python/3.6.0 java bcftools

path2Intersection=$1
normal=$2
panel=$12
csv=$4
alt_cut=$5
tot_cut=$6
vaf_cut=$7
maf_cut=$8
reference=$9
path2database=$10
panelfilter=$11
bam=$3

cd $path2Intersection

if [ -z "$1" ]; then
    err "No path to intersection file provided."
    exit 1
fi

if [ -z "$2" ]; then
    err "No path to germline VCF provided."
    exit 1
fi

normalname="null"

for path in ${path2Intersection}/*; do
    [ -d $path ] || continue

    if [ $panel -eq "False"]; then
        normalname=$(grep "$dirname" ${csv} | cut -d',' -f2 | cut -d'.' -f1)
    else
        normalname=$normal
    fi
    sbatch RunFilter.sh $path2Intersection $normal $csv $alt_cut $tot_cut $vaf_cut $maf_cut $reference $path2database $bam $panelfilter
done

