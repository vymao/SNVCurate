#!/bin/bash

#SBATCH -c 4                               # Request core count
#SBATCH -N 1                               # Request one node (if you request more than one core with -c, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 2-00:00                        # Runtime in D-HH:MM format
#SBATCH -p park                        # Partition to run in
#SBATCH --account=park_contrib
#SBATCH --mem=20000                          # Memory total in MB (for all cores)
#SBATCH -o /n/data1/hms/dbmi/park/victor/other/tests/ReadGroup-log              # File to which STDOUT will be written, including job ID
#SBATCH -e /n/data1/hms/dbmi/park/victor/other/tests/ReadGroup-err           # File to which STDERR will be written, including job ID
#SBATCH --mail-type=FAIL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=victor_mao@hms.harvard.edu   # Email to which notifications will be sent


module load gcc/6.2.0 python/3.6.0 samtools/1.3.1 bwa/0.7.15 java 

path2picard="/home/mk446/BiO/Install/picard-tools-2.5.0/picard.jar"

path=$1
dirname="$(basename "${path}")" 

cd $path

for file in ${path}/*Combined.vcf; do

	bname=$(basename $file)
	base=$(echo $bname | cut -d'.' -f1)

	echo "Splitting..."

  	for sample in `bcftools view -h $file | grep "^#CHROM" | cut -f10-`; do
    	bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
 	done

 	echo 'Indexing...'

	for sample in ${path}/${base}*.gz; do tabix -p vcf -f $sample; done

	echo 'Renaming Sample...'	

	for sample in ${path}/${base}*.gz; do 
	      java -jar $path2picard RenameSampleInVcf \
	      INPUT=$sample\
	      OUTPUT=${sample}.grouped.vcf \
	      NEW_SAMPLE_NAME="TUMOR-Merged"
	done

	echo 'Reindexing...'


	for sample in ${path}/${base}*.gz.grouped.vcf; do 
		bgzip -c $sample > ${sample}.gz
	done


	for sample in ${path}/${base}*.gz.grouped.vcf; do 
		tabix -p vcf -f $sample.gz
	done

	bcftools concat ${base}*.grouped.vcf -o ${base}.Combined.RENAMED.vcf

	rm *.grouped*
	rm *L00*
done
