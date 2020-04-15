#!/bin/bash

#SBATCH -c 1                              # Request core count
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-06:00:00                        # Runtime in D-HH:MM format
#SBATCH -p short                 # Partition to run in
# --account=park_contrib
#SBATCH --mem=1000                          # Memory total in MB (for all cores)
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL


#module load gcc/6.2.0 python/3.6.0 java bcftools

path2SNVCurate=$(dirname "$(readlink -f "$0")")

path2Intersection=$1
normal=$2
matchedNormal=$(echo "$3" | awk '{print tolower($0)}')
csv=$4
alt_cut=$5
tot_cut=$6
vaf_cut=$7
maf_cut=$8
reference=$9
path2database=${10}
bam=${11}
filterwithPanel=$(echo "${12}" | awk '{print tolower($0)}')
panel=${13}


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
    dirname=$(basename $path)

    if [ $matchedNormal == "true" ]; then
        normalname=$(grep "$dirname" ${csv} | cut -d',' -f2 | cut -d'.' -f1)
	base=$(basename $normal)
        dir=$(dirname $normal)
        normalname=${dir}/${base}/${normalname}/${normalname}.vcf
    else
        normalname=$normal
    fi

    if [ $filterwithPanel == "false" ]; then
        sh ${path2SNVCurate}/RunFilter.sh ${path}/intersection_files $normalname $csv $alt_cut $tot_cut $vaf_cut $maf_cut $bam $reference $path2database 
    else
        sh ${path2SNVCurate}/RunFilter.sh ${path}/intersection_files $normalname $csv $alt_cut $tot_cut $vaf_cut $maf_cut $bam $reference $path2database $panel
    fi
done


