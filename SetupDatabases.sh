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

path2database=$1
reference=$2
annovar_database=$3
cd $path2database

ln -s ${annovar_database}/${reference}_refGene.txt ${reference}_refGene.txt
ln -s ${annovar_database}/${reference}_AFR.sites.2015_08.txt ${reference}_AFR.sites.2015_08.txt
ln -s ${annovar_database}/${reference}_ALL.sites.2015_08.txt ${reference}_ALL.sites.2015_08.txt
ln -s ${annovar_database}/${reference}_AMR.sites.2015_08.txt ${reference}_AMR.sites.2015_08.txt
ln -s ${annovar_database}/${reference}_EAS.sites.2015_08.txt ${reference}_EAS.sites.2015_08.txt
ln -s ${annovar_database}/${reference}_EUR.sites.2015_08.txt ${reference}_EUR.sites.2015_08.txt
ln -s ${annovar_database}/${reference}_SAS.sites.2015_08.txt ${reference}_SAS.sites.2015_08.txt
ln -s ${annovar_database}/${reference}_exac03.txt ${reference}_exac03.txt
ln -s ${annovar_database}/${reference}_esp6500siv2_all.txt ${reference}_esp6500siv2_all.txt
ln -s ${annovar_database}/simpleRepeat.bed simpleRepeat.bed
ln -s ${annovar_database}/hg19_rmsk.bed hg19_rmsk.bed
ln -s ${annovar_database}/all_repeats.b37.bed all_repeats.b37.bed
ln -s ${annovar_database}/20141020.strict_mask.whole_genome.bed 20141020.strict_mask.whole_genome.bed
ln -s ${annovar_database}/all.repeatmasker.b37.bed all.repeatmasker.b37.bed
ln -s ${annovar_database}/${reference}_refGeneMrna.fa ${reference}_refGeneMrna.fa
