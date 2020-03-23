# SNVCurate
A filtering pipeline for producing a callset of somatic mutations. This is used only for running on Orchestra, the Harvard Medical School cluster. 

**Note: This pipeline presumes that the BAM files are formatted properly (demultiplexed and indexed correctly). If they are not, please see BAMCurate first before running this pipeline.**

## Dependencies: 
These are the tools/software required (loaded under `module load [dependency]`), which should all be available as modules on O2: 
1. gcc/6.2.0 (the General C Compiler)
2. java (more specifically, the JVM)
3. python/3.6.0
4. bcftools
5. samtools

## Information about relevant scripts: 
1. `Mutect2_read.py`: Wrapper script to run the GATK MuTect2 pipeline for somatic mutation calling. 
```
usage: FastqToSam.py **-tumor [input directory] -out [output directory]**
                     [-n NUM_CORES] [-t RUNTIME] [-p QUEUE]
                     [--mem_per_cpu MEM_PER_CPU] [--mail_type MAIL_TYPE]
                     [--mail_user MAIL_USER] [-picard PICARD_PATH]
                     [-library LIBRARY_NAME]
 ```
                     
    -     parser.add_argument('-tumor', '--input_tumor_path', help='path to input tumor file')
    parser.add_argument('-normal', '--input_normal_path', help='path to normal file')
    parser.add_argument('-pon', '--panel', default='nopath', help='path to panel of norms')
    parser.add_argument('-csv', help='csv containing matched tumor/normal pairs')
    parser.add_argument('-out', '--output_directory', default='./', help='directory to which the output directory "/.Mutect2/" will be written to')
    parser.add_argument('-n', '--num_cores', default='2', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='3-0:0:0', help='slurm job submission option')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='12G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='ALL', help='slurm job submission option')
    parser.add_argument('--mail_user', default='victor_mao@hms.harvard.edu', help='slurm job submission option')
    parser.add_argument('-gatk', '--gatk_path', default='/home/mk446/BiO/Install/GATK3.5_160425_g7a7b7cd/GenomeAnalysisTK.jar', help='path to software')
    # parser.add_argument('-gatk', '--gatk_path', default='/n/data1/hms/dbmi/park/lawrence/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar', help='path to software')
    parser.add_argument('-gatk4', '--gatk4_path', default='/n/data1/hms/dbmi/park/alon/software/gatk/gatk-4.0.3.0/gatk', help='path to software')
    parser.add_argument('-reference', '--reference_path', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta', help='path to reference_path file')
    # parser.add_argument('-reference', '--reference_path', default='/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38/Homo_sapiens_assembly38.fasta', help='path to reference_path 
    parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf', help='path to dbsnp file')
    # parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_hg38_common_all_20160601.vcf', help='path to dbsnp file')
    parser.add_argument('-cosmic', '--cosmic_path', default='/home/mk446/BiO/Install/CosmicCodingMuts_v72.sorted.vcf', help='path to cosmic file' )
    parser.add_argument('-scatter', '--scatter_size', default='200')
    parser.add_argument('-r1', default=1, help='Lower range bound of indices of BAMs to run')
    parser.add_argument('-r2', default=100000, help='Upper range bound of indices of BAMs to run')


## Running MuTect2 (edit parameters as necessary): 
1. `python3 /path/to/Mutect2_read.py -in_dir [input_directory] -out_dir [output_directory] [options]`
    
