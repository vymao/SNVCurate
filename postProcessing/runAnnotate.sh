#!/bin/bash

#SBATCH -c 1                              # Request core count
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-06:00:00                        # Runtime in D-HH:MM format
#SBATCH -p park                 # Partition to run in
#SBATCH --account=park_contrib
#SBATCH --mem=1000                          # Memory total in MB (for all cores)
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL



path=$1
annovarscript=$2
path2database=$3
out=$4
path2Mutect=$5
reference=$6
path2normal=$7
csv=$8
path2SNVCurate=$9


dirname="$(basename "${path}")"
cd ${path}/annotation_files

rm -f *levels
rm -f *LABELED

for file in *somatic_variants_filtered_2.vcf; do 
    outname=${dirname}.somatic_variants_filtered_2
    
    if [ $reference == "hg19" ]; then
        ${annovarscript} $file ${path2database} -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish
    else
        ${annovarscript} $file ${path2database} -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp30a' -operation 'g,f,f' -nastring . -vcfinput -polish
    fi
done

for file in *germline_variants_filtered.vcf; do 
    outname=${dirname}.germline_variants_filtered
    if [ $reference == "hg19" ]; then
        ${annovarscript} $file ${path2database} -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nas$
    else 
        ${annovarscript} $file ${path2database} -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp30a' -operation 'g,f,f' -nas$
    fi
done   

for file in *M2_Risk_variants_filtered.vcf; do 
    outname=${dirname}.M2_Risk_variants_filtered

    if [ $reference == "hg19" ]; then
        ${annovarscript} $file ${path2database} -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nas$
    else 
        ${annovarscript} $file ${path2database} -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp30a' -operation 'g,f,f' -nas$
    fi
done 

if [ -f ${dirname}.PoN_filtered.vcf ]; then 
    outname=${dirname}.PoN_filtered
    
    if [ $reference == "hg19" ]; then
        ${annovarscript} ${dirname}.PoN_filtered.vcf ${path2database} -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish
    else
        ${annovarscript} ${dirname}.PoN_filtered.vcf ${path2database} -out $outname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp30a' -operation 'g,f,f' -nastring . -vcfinput -polish  
    fi
fi


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

if [ -f ${dirname}.PoN_filtered.vcf ]; then 
    for file in *PoN_filtered.*txt; do 
        python3 ${path2SNVCurate}/Label_Source.py -source pon -out ${path}/annotation_files -in $file
    done 

    for file in *PoN_filtered.*txt.LABELED; do 
        python3 ${path2SNVCurate}/Add_Read_Info.py -in_file $file -vcf_path ${dirname}.PoN_filtered.vcf
    done

fi




if [ $path2normal != "False" ]; then
    cd $path2normal

    normalname=$(grep "$dirname" ${csv} | cut -d',' -f2 | cut -d'.' -f1)

    cd ${normalname}*
    
    if [ $reference == "hg19" ]; then
        ${annovarscript} ${path2normal}/${normalname}/${normalname}.vcf ${path2database} -out $normalname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish 
    else
        ${annovarscript} ${path2normal}/${normalname}/${normalname}.vcf ${path2database} -out $normalname -buildver $reference -remove -protocol 'refGene,clinvar_20190305,dbnsfp30a' -operation 'g,f,f' -nastring . -vcfinput -polish
    fi  
     
    if [ ! -f ${normalname}*.LABELED ]; then
        python3 ${path2SNVCurate}/Label_Source.py -source HaplotypeCaller -out ${path2normal}/${normalname}* -in ${normalname}.*txt
    fi
    
    if [ ! -f ${normalname}*.LABELED.levels ]; then
        python3 ${path2SNVCurate}/Add_Read_Info.py -in_file ${normalname}*.LABELED -vcf_path ${path2normal}/${normalname}.vcf -hap True
    fi
   
fi

cd $path

Mutect_Germline_anno_file='*M2_Risk_variants_filtered.*txt.LABELED.levels'
Filtered_file='*germline_variants_filtered.*txt.LABELED.levels'
Haplo_file='.LABELED.levels'
Somatic_file='*somatic_variants_filtered_2.*txt*.levels'
bad_file='ANNO.bad_somatic_quality.*vcf'
pon_file='*PoN_filtered.*.levels'
#Mutect_Germline_risk_file="$dirname.germline_variants_filtered.vcf.TUMOR.avinput.hg19_multianno.csv"

if [ $path2normal != "False" ]; then
    normalname=$(grep "$dirname" ${csv} | cut -d',' -f2 | cut -d'.' -f1)
    if [ ! -f ${out}/${dirname}.germline_combined.csv ]; then
        for file in ${path2normal}/${normalname}/${normalname}*${Haplo_file}; do
            cp $file "${out}/${dirname}.germline_combined.txt"
        done

        for file in ${path}/annotation_files/*${Mutect_Germline_anno_file}; do
            #cp $file "${Output_path}/${dirname}.germline_combined.csv"
            tail -n +2 $file >> "${out}/${dirname}.germline_combined.txt"
        done
        

        for file in ${path}/annotation_files/*${Filtered_file}; do
            tail -n +2 $file >> "${out}/${dirname}.germline_combined.txt"
        done
    fi

else
    normalname=$(grep "$dirname" ${csv} | cut -d',' -f2 | cut -d'.' -f1)
    if [ ! -f ${out}/${dirname}.germline_combined.csv ]; then
        for file in ${path}/annotation_files/*${Mutect_Germline_anno_file}; do
            cp $file "${out}/${dirname}.germline_combined.txt"
       done
        

        for file in ${path}/annotation_files/*${Filtered_file}; do
            tail -n +2 $file >> "${out}/${dirname}.germline_combined.txt"
        done
    fi
fi

if [ -f ${path}/annotation_files/${dirname}.PoN_filtered.vcf ]; then 
    for file in ${path}/annotation_files/*${pon_file}; do
        tail -n +2 $file >> "${out}/${dirname}.germline_combined.txt"
    done
fi

if [ ! -f ${out}/${dirname}.somatic.csv ]; then
    for file in ${path}/annotation_files/*${Somatic_file}; do
        cp $file "${out}/${dirname}.somatic.txt"
    done
fi

if [ ! -f ${out}/${dirname}.bad_somatic_reads.vcf ]; then
    for file in ${path}/cut_filtering/*${bad_file}; do
        cp $file "${out}/${dirname}.bad_somatic_reads.vcf"
    done
fi

