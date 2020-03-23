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


path2run='/n/data1/hms/dbmi/park/victor/scripts/other/Filter_Mutect_Germlines_txt.py'
csv='/n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned/PancSeq_WES_cohort.csv'
path2directory='/n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel/panel_of_normals'


main=$1
path2Mutect=$2
cd $main 
dirname="$(basename "${main}")"

rm *M2_RISK* 
rm *PASS* 
rm *hg19_multianno* 
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
intersected_file=${dirname}.INTERSECTION.MUTECT.vcf
mv 0003.vcf $intersected_file


if [ ! -f ${intersected_file}*txt ]; then
    for file in $intersected_file; do 
        /home/mk446/bin/annovar/table_annovar.pl $file '/home/mk446/bin/annovar/humandb/' -buildver 'hg19' -out $file -remove -protocol 'refGene,avsnp142,exac03,gnomad_genome,1000g2015aug_all' -operation 'g,f,f,f,f' -nastring . -vcfinput -polish
    done
fi


#rm *M2_RISK*
#rm *filtered*
#rm *quality*


for file in *INTERSECTION*.txt; do 
    #echo $file
    mv $file ${file}.csv
done


normalname=$(grep "$dirname" ${csv} | cut -d',' -f4 | cut -d'.' -f1)
#[ -d "${main}" ] || continue 
#echo $dirname
#normal='/n/data1/hms/dbmi/park/victor/references/TCGA_1000_PON.hg19.REORDERED.vcf'
for file in ${main}/*csv; do 
    if [ ! -f ${main}/${dirname}.PASS.vcf.hg19_multianno.txt.ANNO.somatic_variants_filtered.2_5_0.01.vcf ]; then
       #python3 $path2run -input_path ${file} -output_path $main -vcf_path ${main}/${dirname}.PASS.vcf -file_type anno -cut 0.01 -hap $normal -alt_cut 2 -tot_cut 5 -vaf_cut 0.01
       python3 $path2run -input_path ${file} -output_path $main -vcf_path ${main}/${dirname}.INTERSECTION.MUTECT.vcf -file_type anno -cut 0.01 -hap ${path2directory} -alt_cut 2 -tot_cut 5 -vaf_cut 0.01
    fi
done


for file in *ANNO.somatic_variants_filtered.*vcf; do 
    outname=${dirname}.02_05_001.PASS.ANNO.somatic_variants_filtered 
    /home/mk446/bin/annovar/table_annovar.pl $file '/home/mk446/bin/annovar/humandb/' -out $outname -buildver 'hg19' -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish
    sed -i '/#/d' *PASS.ANNO.somatic_variants_filtered*.txt   
done

for file in *ANNO.germline_variants_filtered.*vcf; do 
    outname=${dirname}.02_05_001.PASS.ANNO.germline_variants_filtered
    /home/mk446/bin/annovar/table_annovar.pl $file '/home/mk446/bin/annovar/humandb/' -out $outname -buildver 'hg19' -remove -protocol 'refGene,clinvar_20190305,dbnsfp33a' -operation 'g,f,f' -nastring . -vcfinput -polish
    sed -i '/#/d' *PASS.ANNO.germline_variants_filtered*.txt 
done   


for file in *PASS.ANNO.somatic_variants_filtered.*txt; do 
    python3 /n/data1/hms/dbmi/park/victor/scripts/other/Filter_Mutect_Germlines_Add_Read_Info.py -in_file $file -vcf_path ${main}/${dirname}.INTERSECTION.MUTECT.vcf
done


for file in *PASS.ANNO.germline_variants_filtered.*txt; do 
    python3 /n/data1/hms/dbmi/park/victor/scripts/other/FMG_Label.py -source filtering -out $main -in $file
done 

for file in *PASS.ANNO.germline_variants_filtered.*txt.LABELED; do 
    python3 /n/data1/hms/dbmi/park/victor/scripts/other/Filter_Mutect_Germlines_Add_Read_Info.py -in_file $file -vcf_path ${main}/${dirname}.INTERSECTION.MUTECT.vcf
done 




