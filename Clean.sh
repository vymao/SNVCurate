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


module load gcc/6.2.0 python/3.6.0 java bcftools vcftools


path2MuSE=$3
path2Mutect=$2
path2HaplotypeCaller=$(echo "$1" | awk '{print tolower($0)}')
#normal=$(echo "$2" | awk '{print tolower($0)}')


cd $path2Mutect

for file in ${path2Mutect}/*.vcf; do 
    base=$(basename $file)
    new_dir=$(echo $base | cut -d'.' -f1)
    mkdir -p $new_dir
done

for path in ${path2Mutect}/*; do 
     [ -d $path ] || continue
     sample=$(basename $path)
     mv ${sample}.* $sample
done



if ! [ -z "$path2MuSE" ]; then
     cd $path2MuSE
     for file in ${path2MuSE}/*.vcf; do
	  base=$(basename $file) 
          new_dir=$(echo $base | cut -d'.' -f1)
          mkdir -p $new_dir
     done

     for path in ${path2MuSE}/*; do 
     	[ -d $path ] || continue 
     	sample=$(basename $path)
     	mv ${sample}.* $sample 
     done
fi


if ! [ $path2HaplotypeCaller == "none" ]; then
    cd $path2HaplotypeCaller
    for file in ${path2HaplotypeCaller}/*.vcf; do 
        base=$(basename $file)
	new_dir=$(echo $base | cut -d'.' -f1)
        mkdir -p $new_dir
    done

     for path in ${path2HaplotypeCaller}/*; do
        [ -d $path ] || continue
        sample=$(basename $path)
        mv ${sample}.* $sample
     done
  
    for path in ${path2HaplotypeCaller}/*; do 
        [ -d $path ] || continue
        cd $path
        dirname=$(basename $path)
	if [ ! -f ${dirname}.vcf ]; then
	     bcftools concat *.vcf -o ${dirname}.g.vcf
	fi
	bgzip -c ${dirname}.g.vcf > ${dirname}.g.vcf.gz
	tabix -p vcf ${dirname}.g.vcf.gz
    done
fi



