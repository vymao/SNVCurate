#!/bin/bash


pathtoSNVCurate=$0
path2config=$1


email=$(jq -r '.email' $path2config)
output_directory=$(jq -r '.output_directory' $path2config)
scatter=$(jq -r '.scatter' $path2config)
dbSNP=$(jq -r '.dbSNP' $path2config)
gnomAD=$(jq -r '.gnomAD' $path2config)
reference_path=$(jq -r '.reference_path' $path2config)
panel_path=$(jq -r '.panel_path' $path2config)
reference=$(jq -r '.reference' $path2config)
mail_type=$(jq -r '.mail_type' $path2config)
data_type=$(jq -r '.data_type' $path2config)
csv=$(jq -r '.csv' $path2config)
num_cores=$(jq -r '.num_cores' $path2config)
runtime=$(jq -r '.runtime' $path2config)
queue=$(jq -r '.queue' $path2config)
memory=$(jq -r '.memory' $path2config)
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


cleanCode=$(sh ${pathtoSNVCurate}/Clean.sh $HaplotypeCaller_directory $normal $Mutect_directory $MuSE_directory)

if [ $cleanCode -eq 1 ]; then
    echo "Error in cleaning."
    exit 1
fi

intersectCode=$(sh ${pathtoSNVCurate}/Intersect.sh $output_directory $Mutect_directory $MuSE_directory)

if [ $intersectCode -eq 1 ]; then
    echo "Error in Intersecting."
    exit 1
fi

if [ $normal == "False" ]; then
    filterCode=$(sh ${pathtoSNVCurate}/Filter.sh ${output_directory} $germline_path False $path2bam $csv $alt_cut $tot_cut $vaf_cut $maf_cut $reference $database_path $path2bam $filter_with_panel $panel_path)
else
    filterCode=$(sh ${pathtoSNVCurate}/Filter.sh ${output_directory} $germline_path True $path2bam $csv $alt_cut $tot_cut $vaf_cut $maf_cut $reference $database_path $path2bam $filter_with_panel False)
done


if [ $filterCode -eq 1 ]; then
    echo "Error in Filtering."
    exit 1
fi
out=$1
path2Mutect=$2
reference=$3
path2normal=$5
csv=$4

if [ $normal == "True" ]; then
    annotateCode=$(sh ${pathtoSNVCurate}/Annotate.sh $output_directory $Mutect_directory $reference $csv $germline_path)
else
    annotateCode=$(sh ${pathtoSNVCurate}/Annotate.sh $output_directory $Mutect_directory $reference $csv)
done

if [ $annotateCode -eq 1 ]; then
    echo "Error in Annotating."
    exit 1
fi




