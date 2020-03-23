# SNVCurate
A filtering pipeline for producing a callset of somatic mutations. It also provides a method to curate targeted capture regions in the event that this BED file is not available. 
This pipeline is used only for running on Orchestra, the Harvard Medical School cluster. 

**Note: This pipeline presumes that the BAM files are formatted properly (demultiplexed and indexed correctly). If they are not, please see BAMCurate first before running this pipeline. If no matched normal is available, run this pipeline using another normal of the same sequencing type.**

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
usage: FastqToSam.py [-tumor INPUT_DIRECTORY] [-normal NORMAL_DIRECTORY] [-out OUTPUT_DIRECTORY] [-csv TUMOR/NORMAL_CSV] 
                     [-pon PANEL_OF_NORMALS] [-n NUM_CORES] [-t RUNTIME] [-p QUEUE] [--mem_per_cpu MEM_PER_CPU] [--mail_type MAIL_TYPE]
                     [--mail_user MAIL_USER] [-gatk4 PATH_TO_GATK4] [-reference REFERENCE.FASTA] [-dbsnp dbSNP.vcf] [-gnomad GNOMAD.vcf] 
                     [-scatter] [-r1] [-r2]
```
Additional Information/Default parameters: 
    `-pon`: Used when there are no matched normals. 
    `-csv`: A csv file containing information about matched tumor/normal pairs. See `MuTect2_sample.csv` for proper formatting.
    `-n`: Number of cores (default = 2).
    `-t`: Slurm job runtime. Note that this is the runtime per interval job (default = 0-12:0:0).
    `-p`: Slurm queue (default = park).
    `--mem_per_cpu`: Memory per core (default = 10G).
    `--mail_type`: Notification type (default = FAIL). default='ALL', help='slurm job submission option').
    `--mail_user`: Email.
    `-gatk4`: GATK executable path(default = /n/data1/hms/dbmi/park/alon/software/gatk/gatk-4.0.3.0/gatk).
    `-reference`: Reference FASTA file path (default = /home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta).
    `-dbsnp`: dbSNP database VCF path (default = /home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf).
    `-scatter`: Number of interval files (splits calling into genomic intervals to speed computation) (default = 50).
    `-r1`: The lower range index bound for BAMs to submit from the csv file (default = 1). Index 1 is the lowest. 
    `-r2`: The lower range index bound for BAMs to submit from the csv file (default = 100000).

## Running the SNV curating pipeline: 
There are two ways to run this pipeline: 
1. Input a config file `parameters.config` in the **same directory** as SNVCurate. The config file should follow this format: 

2. Run each component separately. See below for details.

## Running components separately: 
1. Run MuTecT2 on the BAM files using `Mutect2_read.py`. 
2. If you desire higher specifity for somatic calls, run MuSE. 
3. Run the script `Intersect.sh` to intersect the two calls. 
4. Run the script `Filter.sh` to filter the intersection.
5. Run the script `Annotate.sh` to finish. 
    
