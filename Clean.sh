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


main=$4
path2Mutect=$2
path2HaplotypeCaller=$1
normal=$(echo "$1" | awk '{print tolower($0)}')

cd $path2Mutect

for file in ${path2Mutect}/*.vcf; do 
    new_dir=$(echo $file | cut -d'.' -f1)
    if [ ! -d "$new_dir" ]; then
        mkdir $new_dir
    fi
done
for file in ${path2Mutect}/*.vcf; do 
    new_dir=$(echo $file | cut -d'.' -f1)
    mv $file $new_dir
done


cd $main
for file in ${main}/*.vcf; do 
    new_dir=$(echo $file | cut -d'.' -f1)
    if [ ! -d "$new_dir" ]; then
        mkdir $new_dir
    fi
done
for file in ${main}/*.vcf; do 
    new_dir=$(echo $file | cut -d'.' -f1)
    mv $file $new_dir
done

if [ $normal -eq == "true" ]; then
    cd $path2HaplotypeCaller
    for file in ${path2HaplotypeCaller}/*.vcf.gz; do 
        new_dir=$(echo $file | cut -d'.' -f1)
        if [ ! -d "$new_dir" ]; then
            mkdir $new_dir
        fi
    done
    for file in ${path2HaplotypeCaller}/*.vcf.gz; do 
        new_dir=$(echo $file | cut -d'.' -f1)
        mv $file $new_dir
    done
    for path in ${path2HaplotypeCaller}/*; do 
        [ -d $path ] || continue
        cd $path
        dirname=$(basename $path)
        cat *.gz ${dirname}.vcf.gz
        tabix -p vcf ${dirname}.vcf.gz
    done
fi

