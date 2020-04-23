#!/bin/bash


pathtoSNVCurate="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
path2config=$1


email=$(jq -r '.email' $path2config)
output_directory=$(jq -r '.output_directory' $path2config)
reference_path=$(jq -r '.reference_path' $path2config)
panel_path=$(jq -r '.panel_path' $path2config)
reference=$(jq -r '.reference' $path2config)
data_type=$(jq -r '.data_type' $path2config)
csv=$(jq -r '.csv' $path2config)
database_path=$(jq -r '.database_path' $path2config)
Mutect_directory=$(jq -r '.Mutect_directory' $path2config)
MuSE_directory=$(jq -r '.MuSE_path' $path2config)
MuSE=$(jq -r '.MuSE' $path2config)
germline_path=$(jq -r '.germline_path' $path2config)
normal=$(jq -r '.normal' $path2config)
path2bam=$(jq -r '.BAM_path' $path2config)
alt_cut=$(jq -r '.alt_cut' $path2config)
tot_cut=$(jq -r '.tot_cut' $path2config)
vaf_cut=$(jq -r '.vaf_cut' $path2config)
maf_cut=$(jq -r '.maf_cut' $path2config)
filter_with_panel=$(jq -r '.filter_with_panel' $path2config)


set -e 

if [ ${MuSE} == "True" ]; then
    intersectCode=$(sh ${pathtoSNVCurate}/Intersect.sh $output_directory $Mutect_directory $MuSE_directory)
else
    intersectCode=$(sh ${pathtoSNVCurate}/Intersect.sh $output_directory $Mutect_directory)
fi

filterCode=$(sh ${pathtoSNVCurate}/Filter.sh ${output_directory} $germline_path $normal $csv $alt_cut $tot_cut $vaf_cut $maf_cut $reference $database_path $path2bam $filter_with_panel $panel_path)


if [ $normal == "True" ]; then
    annotateCode=$(sh ${pathtoSNVCurate}/Annotate.sh $output_directory $Mutect_directory $reference $csv $germline_path)
else
    annotateCode=$(sh ${pathtoSNVCurate}/Annotate.sh $output_directory $Mutect_directory $reference $csv)
fi




