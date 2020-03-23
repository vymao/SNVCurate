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


main=$1
path2Mutect=$2

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



cd $path2Mutect
for path in ${path2Mutect}/*; do 
    [ -d "${path}" ] || continue 
    cd $path

    dirname="$(basename "${path}")"

    rm *PASS* 
    rm *INTERSECTION* *TIER* 000* *MUTECT*

    if [ ! -f ${dirname}.PASS.vcf ]; then
        for file in ${dirname}.vcf; do 
            grep "PASS\|#" $file > ${dirname}.PASS.vcf
        done
    fi
done

cd $main 
for path in ${main}/*; do 
    [ -d "${path}" ] || continue 
    cd $path

    dirname="$(basename "${path}")"

    rm *PASS* 
    rm *INTERSECTION* *TIER* 000* *MUTECT*

    if [ ! -f ${dirname}.PASS.vcf ]; then
        for file in ${dirname}.vcf; do 
            grep "PASS\|#" $file > ${dirname}.PASS.vcf
        done
    fi

    if [ ! -f ${dirname}.TIER.vcf ]; then
        for file in ${dirname}.vcf; do 
            grep "Tier\|#" $file > ${dirname}.TIER.vcf
        done
    fi
done

bcftools sort ${dirname}.PASS.vcf > ${dirname}.PASS.sorted.vcf
cp ${dirname}.PASS.sorted.vcf ${dirname}.PASS.vcf

bgzip -c ${dirname}.PASS.vcf > ${dirname}.PASS.vcf.gz
tabix -p vcf ${dirname}.PASS.vcf.gz

mutectFile=${path2Mutect}/${dirname}/${dirname}.PASS.vcf
bcftools sort $mutectFile > ${dirname}.MUTECT_SORTED.vcf
mutectFile=${dirname}.MUTECT_SORTED.vcf
bgzip -c $mutectFile > ${mutectFile}.gz
tabix -p vcf ${mutectFile}.gz

bcftools isec -p $PWD -Oz ${dirname}.PASS.vcf.gz ${mutectFile}.gz
gunzip 0003.vcf.gz
intersected_file=${dirname}.INTERSECTION.vcf
mv 0003.vcf $intersected_file
