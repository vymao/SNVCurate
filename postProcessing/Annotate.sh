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

rm *csv
rm *vcf

for path in ${out}/*; do 
    [ -d $path ] || continue
    dirname="$(basename "${path}")"
    cd ${path}/annotation_files

    rm *levels
    rm *LABELED

    for file in *somatic_variants_filtered_2.vcf; do 
        outname=${dirname}.somatic_variants_filtered_2
        /home/mk446/bin/annovar/table_annovar.pl $file '/home/mk446/bin/annovar/humandb/' -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish
        #sed -i '/#/d' *PASS.ANNO.somatic_variants_filtered*.txt   
    done

    for file in *germline_variants_filtered.vcf; do 
        outname=${dirname}.germline_variants_filtered
        /home/mk446/bin/annovar/table_annovar.pl $file '/home/mk446/bin/annovar/humandb/' -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish
        #sed -i '/#/d' *PASS.ANNO.germline_variants_filtered*.txt 
    done   

    for file in *M2_Risk_variants_filtered.vcf; do 
        outname=${dirname}.M2_Risk_variants_filtered
        /home/mk446/bin/annovar/table_annovar.pl $file '/home/mk446/bin/annovar/humandb/' -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish
        #sed -i '/#/d' ${dirname}.M2_RISK.germline_variants_filtered*.txt 
    done 


    mutect=${out}/${dirname}/intersection_files


    for file in *somatic_variants_filtered_2.*txt; do 
        python3 ${path2SNVCurate}/Add_Read_Info.py -in_file $file -vcf_path ${mutect}/*PASS_MuTecT.vcf
    done

    for file in *germline_variants_filtered.*txt; do 
        python3 ${path2SNVCurate}/Label_Source.py -source filtering -out ${path}/annotation_files -in $file
    done 

    for file in *germline_variants_filtered.*txt.LABELED; do 
        python3 ${path2SNVCurate}/Add_Read_Info.py -in_file $file -vcf_path ${mutect}/*PASS_MuTecT.vcf
    done 

    for file in *M2_Risk_variants_filtered.*txt; do 
       python3 ${path2SNVCurate}/Label_Source.py -source Mutect -out ${path}/annotation_files -in $file
    done 

    for file in *M2_Risk_variants_filtered.*txt.LABELED; do 
        python3 ${path2SNVCurate}/Add_Read_Info.py -in_file $file -vcf_path ${path2Mutect}/${dirname}/${dirname}.vcf
    done 



    if ! [ -z "$path2normal" ]; then
        cd $path2normal

        dirname=$(grep "$dirname" ${csv} | cut -d',' -f2 | cut -d'.' -f1)

        cd ${dirname}*

        /home/mk446/bin/annovar/table_annovar.pl ${path2normal}/${dirname}/${dirname}.vcf '/home/mk446/bin/annovar/humandb/' -out $dirname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish

         
        if [ ! -f ${dirname}*.LABELED ]; then
            python3 ${path2SNVCurate}/Label_Source.py -source HaplotypeCaller -out ${path2normal}/${dirname}* -in ${dirname}.*txt
        fi
        
        if [ ! -f ${dirname}*.LABELED.levels ]; then
            python3 ${path2SNVCurate}/Add_Read_Info.py -in_file ${dirname}*.LABELED -vcf_path ${path2normal}/${dirname}.vcf -hap True
        fi
       
    fi


done



for path in ${out}/*; do
    [ -d "${path}" ] || continue # if not a directory, skip'
    dirname="$(basename "${path}")"
    cd $path

    Mutect_Germline_anno_file='*M2_Risk_variants_filtered.*txt.LABELED.levels'
    Filtered_file='*germline_variants_filtered.*txt.LABELED.levels'
    Haplo_file='.LABELED.levels'
    Somatic_file='*somatic_variants_filtered_2.*txt*.levels'
    bad_file='ANNO.bad_somatic_quality.*vcf'
    #Mutect_Germline_risk_file="$dirname.germline_variants_filtered.vcf.TUMOR.avinput.hg19_multianno.csv"

    if ! [ -z "$path2normal" ]; then
        normalname=$(grep "$dirname" ${csv} | cut -d',' -f2 | cut -d'.' -f1)
        if [ ! -f ${out}/${dirname}.germline_combined.csv ]; then
            for file in ${path2normal}/${normalname}/${normalname}*${Haplo_file}; do
                cp $file "${out}/${dirname}.germline_combined.csv"
            done

            for file in ${path}/annotation_files/*${Mutect_Germline_anno_file}; do
                #cp $file "${Output_path}/${dirname}.germline_combined.csv"
                tail -n +2 $file >> "${out}/${dirname}.germline_combined.csv"
            done
            

            for file in ${path}/annotation_files/*${Filtered_file}; do
                tail -n +2 $file >> "${out}/${dirname}.germline_combined.csv"
            done
        fi

    else
        normalname=$(grep "$dirname" ${csv} | cut -d',' -f2 | cut -d'.' -f1)
        if [ ! -f ${out}/${dirname}.germline_combined.csv ]; then
            for file in ${path}/annotation_files/*${Mutect_Germline_anno_file}; do
                cp $file "${out}/${dirname}.germline_combined.csv"
           done
            

            for file in ${path}/annotation_files/*${Filtered_file}; do
                tail -n +2 $file >> "${out}/${dirname}.germline_combined.csv"
            done
        fi
    fi


    if [ ! -f ${out}/${dirname}.somatic.csv ]; then
        for file in ${path}/annotation_files/*${Somatic_file}; do
            cp $file "${out}/${dirname}.somatic.csv"
        done
    fi

    if [ ! -f ${out}/${dirname}.bad_somatic_reads.vcf ]; then
        for file in ${path}/cut_filtering/*${bad_file}; do
            cp $file "${out}/${dirname}.bad_somatic_reads.vcf"
        done
    fi

    var=$((var+1))   

done


