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


path2Intersection=$1
normal=$2
csv=$3 
alt_cut=$4 
tot_cut=$5 
vaf_cut=$6 
maf_cut=$7

sampledir=$(dirname $path2Intersection)
dirname=$(basename $sampledir)

cd $sampledir
mkdir cut_filtering
mkdir pon_filtering


cd cut_filtering
if [ ! -f ${path2Intersection}*txt ]; then
    /home/mk446/bin/annovar/table_annovar.pl ${path2Intersection}/${dirname}.INTERSECTION.vcf '/home/mk446/bin/annovar/humandb/' -buildver 'hg19' -out $dirname -remove -protocol 'refGene,avsnp142,exac03,gnomad_genome,1000g2015aug_all' -operation 'g,f,f,f,f' -nastring . -vcfinput -polish
fi

for file in *.txt; do 
    mv $file ${file}.csv
done

for file in *csv; do 
    if [ ! -f cut_filtering/${dirname}.PASS.vcf.hg19_multianno.txt.ANNO.somatic_variants_filtered.*.vcf ]; then
       python3 Filter_Mutect_Germlines_txt.py -input_path ${file} -output_path ${sampledir}/cut_filtering -vcf_path ${path2Intersection}/${dirname}.INTERSECTION.vcf -file_type anno -cut $maf_cut -hap $normal -alt_cut $alt_cut -tot_cut $tot_cut -vaf_cut $vaf_cut
    fi
done

mv ${dirname}.PASS.vcf.hg19_multianno.txt.ANNO.somatic_variants_filtered.*.vcf ${dirname}.somatic_variants_filtered_1.vcf







