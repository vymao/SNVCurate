#!/bin/bash

#SBATCH -c 1                              # Request core count
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-06:00:00                        # Runtime in D-HH:MM format
#SBATCH -p park                 # Partition to run in
#SBATCH --account=park_contrib
#SBATCH --mem=10000                          # Memory total in MB (for all cores)
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL


#module load gcc/6.2.0 python/3.6.0 java bcftools

path2SNVCurate=${12}
echo $0

path2Intersection=$1
normal=$2
csv=$3 
alt_cut=$4 
tot_cut=$5 
vaf_cut=$6 
maf_cut=$7
panelfilter=${11}
reference=$9
path2database=${10}
bam=$8
path2AnnovarScript=${13}


sampledir=$(dirname $path2Intersection)
dirname=$(basename $sampledir)

cd $sampledir
rm -Rf cut_filtering
mkdir -p cut_filtering
#mkdir pon_filtering


cd cut_filtering 

if [ ! -f ${path2Intersection}*txt ]; then
    if [ $reference == "hg19" ]; then
        ${path2AnnovarScript} ${path2Intersection}/${dirname}.INTERSECTION.vcf ${path2database} -buildver ${reference} -out $dirname -remove -protocol 'refGene,exac03,g$
    else
        ${path2AnnovarScript} ${path2Intersection}/${dirname}.INTERSECTION.vcf ${path2database} -buildver ${reference} -out $dirname -remove -protocol 'refGene,exac03,gnomad_genome,gnomad_exome,1000g2015aug_all' -operation 'g,f,f,f,f' -nastring . -vcfinput -polish
    fi
fi

for file in *.txt; do 
    mv $file ${file}.csv
done

set -e

for file in *csv; do 
    if [ ! -f cut_filtering/${dirname}.PASS.vcf.hg19_multianno.txt.ANNO.somatic_variants_filtered.*.vcf ]; then
       echo "Running command python3 ${path2SNVCurate}/Filter_Mutect_Germlines_txt.py -input_path ${file} -output_path ${sampledir}/cut_filtering -vcf_path ${path2Intersection}/${dirname}.INTERSECTION.vcf -file_type anno -cut $maf_cut -hap $normal -alt_cut $alt_cut -tot_cut $tot_cut -vaf_cut $vaf_cut"
       python3 ${path2SNVCurate}/Filter_Mutect_Germlines_txt.py -input_path ${file} -output_path ${sampledir}/cut_filtering -vcf_path ${path2Intersection}/${dirname}.INTERSECTION.vcf -file_type anno -cut $maf_cut -hap $normal -alt_cut $alt_cut -tot_cut $tot_cut -vaf_cut $vaf_cut
    fi
done

mv ${dirname}.*.ANNO.somatic_variants_filtered* ${dirname}.somatic_variants_filtered_1.vcf
mv ${dirname}.*.ANNO.germline_variants_filtered* ${dirname}.germline_variants_filtered.vcf

if [ ${panelfilter} != "False" ]; then
    cd ${sampledir}/cut_filtering
    echo "Running with parameters:"
    echo "alt_cut: ${alt_cut}"
    echo "total_cut: ${tot_cut}"
    echo "VAF_cut: ${vaf_cut}"
    echo "MAF_cut: ${maf_cut}"

    python3 ${path2SNVCurate}/PoN_filter.py -somatic_vcf ${sampledir}/cut_filtering/${dirname}.somatic_variants_filtered_1.vcf -normal_vcf $normal -annovar $path2database -reference $reference -bam $bam -pon $panelfilter
    
    mv ${sampledir}/cut_filtering/${dirname}.somatic_variants_filtered_1.vcf ${sampledir}/annotation_files
    mv ${dirname}_Final_Callset.vcf ${sampledir}/annotation_files/${dirname}.somatic_variants_filtered_2.vcf
    mv ${dirname}_Filtered_file.vcf ${sampledir}/annotation_files/${dirname}.PoN_filtered.vcf
    mv ${dirname}.germline_variants_filtered.vcf ${sampledir}/annotation_files
else 
    cp ${dirname}.somatic_variants_filtered_1.vcf ${sampledir}/annotation_files/${dirname}.somatic_variants_filtered_2.vcf
    mv ${dirname}.somatic_variants_filtered_1.vcf ${sampledir}/annotation_files
    mv ${dirname}.germline_variants_filtered.vcf ${sampledir}/annotation_files
fi




