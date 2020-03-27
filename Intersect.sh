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


main=$3
path2Mutect=$2
out=$1

cd $path2Mutect
for path in ${path2Mutect}/*; do 
    [ -d "${path}" ] || continue 
    cd $path

    dirname="$(basename "${path}")"

    rm *PASS* 
    rm *INTERSECTION* *TIER* 000* *MUTECT*

    bcftools concat *.vcf -o ${dirname}.Combined.vcf

    if [ ! -f ${dirname}.PASS_MuTecT.vcf ]; then
        for file in ${dirname}.Combined.vcf; do 
            grep "PASS\|#" $file > ${dirname}.PASS_MuTecT.vcf
        done
    fi

    for file in ${dirname}/*Combined.vcf; do 
        python3 Filter_Mutect_Germlines_txt.py -input_path $file -output_path ${path} -vcf $file -file_type vcf 
    done

    test_dir=${out}/${dirname}
    mkdir $test_dir
    mkdir ${test_dir}/annotation_files
    mv ${dirname}.M2_RISK.germline_variants_filtered.vcf ${dirname}.M2_Risk_variants_filtered.vcf
    mv ${dirname}.M2_Risk_variants_filtered.vcf ${test_dir}/annotation_files

    mkdir ${test_dir}/intersection_files
    mv ${dirname}.PASS_MuTecT.vcf ${test_dir}/intersection_files
done

if ! [ -z "$main" ]; then
    cd $main 
    for path in ${main}/*; do 
        [ -d "${path}" ] || continue 
        cd $path

        dirname="$(basename "${path}")"

        rm *PASS* 
        rm *INTERSECTION* *TIER* 000* *MUTECT*

        if [ ! -f ${dirname}.PASS_MuSE.vcf ]; then
            for file in ${dirname}.vcf; do 
                grep "PASS\|#" $file > ${dirname}.PASS_MuSE.vcf
            done
        fi

        if [ ! -f ${dirname}.TIER.vcf ]; then
            for file in ${dirname}.vcf; do 
                grep "Tier\|#" $file > ${dirname}.TIER.vcf
            done
        fi

        test_dir=${out}/${dirname}
        mkdir $test_dir
        mkdir ${test_dir}/intersection_files
        mv ${dirname}.PASS_MuSE.vcf $test_dir/intersection_files
    done


    for path in $out/*; do 
        [ -d $path ] || continue
        cd $path
        dirname=$(basename $path)
        cd intersection_files

        bcftools sort ${dirname}.PASS_MuSE.vcf > ${dirname}.PASS_MuSE.sorted.vcf
        cp ${dirname}.PASS_MuSE.sorted.vcf ${dirname}.PASS_MuSE.vcf

        bgzip -c ${dirname}.PASS_MuSE.vcf > ${dirname}.PASS_MuSE.vcf.gz
        tabix -p vcf ${dirname}.PASS_MuSE.vcf.gz

        mutectFile=${path2Mutect}/${dirname}/${dirname}.PASS.vcf
        bcftools sort ${dirname}.PASS_MuTecT.vcf > ${dirname}.PASS_MuTecT.sorted.vcf 
        cp ${dirname}.PASS_MuSE.sorted.vcf ${dirname}.PASS_MuTecT.vcf

        bgzip -c ${dirname}.PASS_MuTecT.vcf > ${dirname}.PASS_MuTecT.vcf.gz
        tabix -p vcf ${dirname}.PASS_MuTecT.vcf.gz

        bcftools isec -p $PWD -Oz ${dirname}.PASS_MuSE.vcf.gz ${dirname}.PASS_MuTecT.vcf.gz
        gunzip 0003.vcf.gz
        intersected_file=${dirname}.INTERSECTION.vcf
        mv 0003.vcf $intersected_file
    done
else
    for path in ${path2Mutect}/*; do 
    [ -d "${path}" ] || continue 
    cd ${path}/intersection_files

    dirname="$(basename "${path}")"
    intersected_file=${dirname}.INTERSECTION.vcf

    mv ${dirname}.PASS_MuTecT.vcf $intersected_file
done
fi







