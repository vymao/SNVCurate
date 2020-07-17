#!/bin/bash

#SBATCH -c 1                              # Request core count
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00:00                        # Runtime in D-HH:MM format
#SBATCH -p medium                 # Partition to run in
#SBATCH --account=park_contrib
#SBATCH --mem=1000                          # Memory total in MB (for all cores)
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL


#module load gcc/6.2.0 python/3.6.0 java bcftools

path2SNVCurate="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

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
annovarscript=${12}
filterwithPanel=$(echo "${13}" | awk '{print tolower($0)}')
panel=${14}


cd $path2Intersection

if [ ! ${matchedNormal} == "false" ]; then
    if [ ! ${matchedNormal} == "true" ]; then
        err "Matched normal Boolean field must be True or False"
        exit 1
    fi
fi

if [ ! ${filterwithPanel} == "false" ]; then
    if [ ! ${filterwithPanel} == "true" ]; then
        err "Filter with masks Boolean field must be True or False"
        exit 1
    fi
fi


normalname="null"

for path in ${path2Intersection}/*; do
    [ -d $path ] || continue
    cd $path
    mkdir -p batch_submissions
    cd batch_submissions

    dirname=$(basename $path)
    ln -s ${path2SNVCurate}/RunFilter.sh ${dirname}_FILTER.sh

    if [ $matchedNormal == "true" ]; then
        normalname=$(grep "$dirname" ${csv} | cut -d',' -f2 | cut -d'.' -f1)
        base=$(basename $normal)
        dir=$(dirname $normal)
        normalname=${dir}/${base}/${normalname}/${normalname}.vcf.gz
    else
        normalname=$normal
    fi

    if [ $filterwithPanel == "false" ]; then
        echo $matchedNormal
        sbatch ${dirname}_FILTER.sh ${path}/intersection_files $normalname $csv $alt_cut $tot_cut $vaf_cut $maf_cut $bam $reference $path2database False $path2SNVCurate $annovarscript $matchedNormal
    else
        echo $matchedNormal
        sbatch ${dirname}_FILTER.sh ${path}/intersection_files $normalname $csv $alt_cut $tot_cut $vaf_cut $maf_cut $bam $reference $path2database $panel $path2SNVCurate $annovarscript $matchedNormal
    fi
done


