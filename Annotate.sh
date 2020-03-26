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

out=$1
path2Mutect=$2
reference=$3
path2normal=$5
csv=$4

cd $out

for path in ${out}/*; do 
    [ -d $path ] || continue
    dirname="$(basename "${path}")"
    cd annotation_files

    for file in *somatic_variants_filtered_2.vcf; do 
        outname=${dirname}.somatic_variants_filtered_1
        /home/mk446/bin/annovar/table_annovar.pl $file '/home/mk446/bin/annovar/humandb/' -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish
        sed -i '/#/d' *PASS.ANNO.somatic_variants_filtered*.txt   
    done

    for file in *germline_variants_filtered.vcf; do 
        outname=${dirname}.02_05_001.PASS.ANNO.germline_variants_filtered
        /home/mk446/bin/annovar/table_annovar.pl $file '/home/mk446/bin/annovar/humandb/' -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish
        sed -i '/#/d' *PASS.ANNO.germline_variants_filtered*.txt 
    done   

    for file in *M2_Risk_variants_filtered.vcf; do 
        outname=${dirname}.M2_RISK.germline_variants_filtered
        /home/mk446/bin/annovar/table_annovar.pl $file '/home/mk446/bin/annovar/humandb/' -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish
        sed -i '/#/d' ${dirname}.M2_RISK.germline_variants_filtered*.txt 
    done 


    mutect=${path2Mutect}/${dirname}


    for file in *somatic_variants_filtered_2.*txt; do 
        python3 Add_Read_Info.py -in_file $file -vcf_path ${mutect}/*PASS.vcf
    done


    for file in *germline_variants_filtered.*txt; do 
        python3 Label_Source.py -source filtering -out $out -in $file
    done 

    for file in *germline_variants_filtered.*txt.LABELED; do 
        python3 Add_Read_Info.py -in_file $file -vcf_path ${mutect}/*PASS.vcf
    done 

    for file in *M2_Risk_variants_filtered.*txt; do 
       python3 Label_Source.py -source Mutect -out $out -in $file
    done 

    for file in *M2_Risk_variants_filtered.*txt.LABELED; do 
        python3 Add_Read_Info.py -in_file $file -vcf_path ${mutect}/*.Combined.vcf
    done 



    if ! [ -z "$path2normal" ]; then
        cd $path2normal

        /home/mk446/bin/annovar/table_annovar.pl ${path2normal}/${dirname}.vcf '/home/mk446/bin/annovar/humandb/' -out $dirname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish

         
        if [ ! -f ${dirname}*.LABELED ]; then
            python3 Label_Source.py -source HaplotypeCaller -out $path2normal -in ${dirname}.*txt
        fi
        
        if [ ! -f ${dirname}*.LABELED.levels ]; then
            python3 Add_Read_Info.py -in_file ${dirname}*.LABELED -vcf_path ${path2normal}/${dirname}.vcf -hap True
        fi
       
    fi


done



for path in $out/*; do
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
            for file in $path2normal/${normalname}*${Haplo_file}; do
                cp $file "${out}/${dirname}.germline_combined.csv"
            done

            for file in ${path}/*${Mutect_Germline_anno_file}; do
                #cp $file "${Output_path}/${dirname}.germline_combined.csv"
                tail -n +2 $file >> "${Output_path}/${dirname}.germline_combined.csv"
            done
            

            for file in ${path}/*${Filtered_file}; do
                tail -n +2 $file >> "${Output_path}/${dirname}.germline_combined.csv"
            done
        fi

    else
        normalname=$(grep "$dirname" ${csv} | cut -d',' -f2 | cut -d'.' -f1)
        if [ ! -f ${out}/${dirname}.germline_combined.csv ]; then
            for file in ${path}/*${Mutect_Germline_anno_file}; do
                cp $file "${Output_path}/${dirname}.germline_combined.csv"
            done
            

            for file in ${path}/*${Filtered_file}; do
                tail -n +2 $file >> "${Output_path}/${dirname}.germline_combined.csv"
            done
        fi
    fi


    if [ ! -f ${out}/${dirname}.somatic.csv ]; then
        for file in ${path}/*${Somatic_file}; do
            cp $file "${out}/${dirname}.somatic.csv"
        done
    fi

    if [ ! -f ${out}/${dirname}.bad_somatic_reads.vcf ]; then
        for file in ${path}/*${bad_file}; do
            cp $file "${out}/${dirname}.bad_somatic_reads.vcf"
        done
    fi

    var=$((var+1))   

done


