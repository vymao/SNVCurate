# SNVCurate
A pipeline for producing a callset of somatic mutations. It also provides a method to curate targeted capture regions in the event that this BED file is not available. 
This pipeline is used only for running on Orchestra, the Harvard Medical School cluster, with the SLURM job scheduler. This pipeline also utilizes GATK 4.1.2.0, the latest version of GATK available on the cluster.

**Note: This pipeline presumes that the BAM files are formatted properly (demultiplexed and indexed correctly). If no matched normal is available, run this pipeline using another normal of the same sequencing type.**

The pipeline is split into two parts: one to generate callsets from Mutect2 (and MuSE + HaplotypeCaller if normals are available), and one to filter those calls. Callset generation is automated on the cluster using the Cromwell execution engine and customized .wdl scripts produced from the datasets inputted. Post-processing is automated using bash wrapper scripts utilizing various software pre-loaded in the environment. It is split this way to prevent overloading of the job partitions, so you can continue submitting jobs as they complete.

If you have some samples with matched normals and some without that you would like to call, you should create separate .csv files for each of these datasets and run this pipeline on each set individually. 

If you would like to run the filtering portion only on existing VCFs, then please see the Wiki for how to do so. 

## Running the SNV curating pipeline (see below for more information on individual scripts): 
### Environment/file setup:
1. Create two files: One csv with tumor/normal matched pairs (leave `N/A` if no normal) and one text file with each normal sample (if applicable). See `Tumor_Normal_sample.csv` and `HaplotypeCaller_sample.txt` for formatting.
2. Set up the environment using Conda, or download and install the packages individually (see .yml file for dependencies). Other modules available on the cluster will be loaded automatically within the relevant script.
3. Run `RenameBAMsample.sh` to reconfigure the BAM sample name to single tumor/normal. You should use this directory as the source of BAM files for later steps.
4. Run `SetupDatabases.sh` to setup links to relevant Annovar databases. 
5. Activate the environment. 

### Calling (under `SNVCurate/calling/`): 
1. Run Mutect2 on the renamed BAM files using `Mutect2_read.py`. 
2. If you desire higher specifity for somatic calls (recommended), run MuSE using `MuSE_read.py`. 
3. If there are matched normals, run HaplotypeCaller on these normals using `HaplotypeCaller_read.py`.

### Post-Processing (under `SNVCurate/postProcessing/`): 
This is best run on an interactive session with 5G of memory. 

1. Run the script `Intersect.sh` to intersect the two calls. 
2. Run the script `Filter.sh` to filter the intersection.
3. Run the script `Annotate.sh` to finish. 

Note that `Filter.sh` and `Annotate.sh` will submit SLURM batch jobs, which you should wait for to finish until moving to the next step. To change the parameters of the SLURM batch jobs (ie. time, queue, node count, etc.), change the parameters of the headers of the scripts `runRenaming.sh`, `/postProcessing/runAnnotate.sh`, and `/postProcessing/RunFilter.sh`. 

## Example:
The files under `SNVCurate/test/` were used for generating the callsets within that folder. The BAM files used for calling are at `/n/data1/hms/dbmi/park/victor/other/pipeline_test/bam_files`, and the steps utilized documented below (the actual file paths should be changed). 

Running an interactive session: 
```
srun --pty -t 0-2:0:0 --mem 5G -p interactive /bin/bash
```

For setup and calling: 
```
conda env create -f environment.yml
conda activate SNVCurate
sh RenameBAMsample.sh /path/to/bam_files /path/to/empty/directory/for/renamed_bams
sh SetupDatabases.sh /path/to/empty/directory/for/databases hg19

python3 Mutect2_read.py -tumor /path/to/bam_files -normal /path/to/bam_files -out /path/to/Mutect2_output_directory -csv /path/to/tumor_normal.csv --mail_user victor_mao@hms.harvard.edu -r1 1 -r2 2 -p park
python3 MuSE_read.py -tumor -tumor /path/to/bam_files -normal /path/to/bam_files -out /path/to/MuSE_output_directory -csv /path/to/tumor_normal.csv --mail_user victor_mao@hms.harvard.edu -data_type WES -r1 1 -r2 2 -p medium
python3 HaplotypeCaller_read.py -csv /path/to/tumor_normal.csv -normal /path/to/bam_files -out /path/to/HaplotypeCaller_output_directory --mail_user victor_mao@hms.harvard.edu -r1 1 -r2 2 -p park -reference_name b37
```

For the filtering:
```
sh Intersect.sh /path/to/filtering_output_directory /path/to/Mutect2_output_directory /path/to/MuSE_output_directory
sh Filter.sh /path/to/filtering_output_directory /path/to/HaplotypeCaller_output_directory True /path/to/tumor_normal.csv 4 10 0.05 0.01 hg19 /path/to/databases /path/to/bam_files True /path/to/panel
sh Annotate.sh /path/to/filtering_output_directory /path/to/Mutect2_output_directory hg19 /path/to/tumor_normal.csv /path/to/HaplotypeCaller_output_directory 
```

## Information about relevant scripts: 

1. `RenameBAMsample.sh`: Bash script to rename `SM` tag in the read groups for each sample according to the sample name. Note: this takes quite a bit of memory and computational capacity so running on an interactive session with sufficient memory is recommended.
```
usage: sh RenameBAMsample.sh [BAM_DIRECTORY] [OUTPUT_DIRECTORY] 
```
- `[BAM_DIRECTORY]`: Directory for BAM files. 
- `[OUTPUT_DIRECTORY]`: Output directory for renamed BAM files. 

2. `SetupDatabases.sh`: Bash script to create links to relevant databases and set up a space for new accessible databases.  
```
usage: sh SetupDatabases.sh [OUTPUT_DIRECTORY] [REFERENCE] [ANNOVAR_DATABASES]
```
- `[OUTPUT_DIRECTORY]`: Output directory for renamed BAM files. 
- `[REFERENCE]`: hg19 or hg38 (for Annovar). 
- `[ANNOVAR_DATABASE]`: Path to Annovar databases. On Orchestra, this is `/home/mk446/bin/annovar/humandb/`. 

3. `Mutect2_read.py`: Wrapper script to run the GATK MuTect2 pipeline for somatic mutation calling. 
```
usage: python3 Mutect2_read.py [-tumor INPUT_DIRECTORY] [-normal NORMAL_DIRECTORY] [-out OUTPUT_DIRECTORY] [-csv TUMOR/NORMAL_CSV] 
                     [-pon PANEL_OF_NORMALS] [-n NUM_CORES] [-t RUNTIME] [-p QUEUE] [--mem_per_cpu MEM_PER_CPU] [--mail_type MAIL_TYPE]
                     [--mail_user MAIL_USER] [-gatk PATH_TO_GATK4] [-reference REFERENCE.FASTA] [-dbsnp dbSNP.vcf] [-gnomad GNOMAD.vcf] 
                     [-scatter COUNT] [-cn CROMWELL_CORES] [-ct CROMWELL_RUNTIME] [-cm CROMWELL_MEMORY] [-cromwell CROMWELL_JAR] 
                     [-interval_list INTERVAL_LIST] [-parallel RUN_IN_PARALLEL] [-r1 START] [-r2 END]
```
- `-pon`: Used when there are no matched normals. 
- `-csv`: A csv file containing information about matched tumor/normal pairs. See `MuTect2_sample.csv` for proper formatting.
- `-n`: Number of cores (default = 1).
- `-t`: Slurm job runtime. Note that this is the runtime per interval job (default = 0-12:0:0).
- `-p`: Slurm queue (default = park).
- `--mem_per_cpu`: Memory per core (default = 10G).
- `--mail_type`: Notification type (default = FAIL)
- `--mail_user`: Email.
- `-gatk`: GATK executable path(default = `/n/data1/hms/dbmi/park/alon/software/gatk/gatk-4.0.3.0/gatk`).
- `-reference`: Reference FASTA file path (default = `/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta`). If different from default, `Mutect_read.py` will produce a text file containing a new list of interval files. 
- `-dbsnp`: dbSNP database VCF path (default = `/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf`).
- `-scatter`: Number of interval files (splits calling into genomic intervals to speed computation) (default = 50). This will determine the number of Cromwell jobs scattered; if different from 50, `Mutect_read.py` will produce a text file containing a new list of interval files. 
- `-interval_list`: Interval list used to scatter Cromwell jobs. The default is set to scatter 50 jobs (default = `/n/data1/hms/dbmi/park/victor/software/MuTecT2_b37_scattered_intervals.txt`) for the b37 reference.
- `-r1`: The lower range index bound for BAMs to submit from the csv file (default = 1). Index 1 is the lowest. 
- `-r2`: The upper range index bound for BAMs to submit from the csv file (default = 100000).
- `-cn`: Number of cores per Cromwell job (default = 1). 
- `-ct`: Runtime per Cromwell job in minutes (default = 1000). 
- `-cm`: Memory per Cromwell job in MB (default = 7000). 
- `-cromwell`: Jar file for Cromwell execution (default = `/n/shared_db/singularity/hmsrc-gatk/cromwell-43.jar`).
- `-parallel`: Parallelize the pipeline for calling in parallel. Only available for matched normal calling (default = `True`).

4. `MuSE_read.py`: Wrapper script to run the MuSE pipeline for somatic mutation calling. 
```
usage: python3 MuSE_read.py [-tumor INPUT_DIRECTORY] [-normal NORMAL_DIRECTORY] [-out OUTPUT_DIRECTORY] [-csv TUMOR/NORMAL_CSV] 
                            [-n NUM_CORES] [-t RUNTIME] [-p QUEUE] [--mem_per_cpu MEM_PER_CPU] [--mail_type MAIL_TYPE] 
                            [--mail_user MAIL_USER] [-reference REFERENCE.FASTA] [-dbsnp dbSNP.vcf.gz] [-data_type WGS/WES] 
                            [-r1 START] [-r2 END] [-cn CROMWELL_CORES] [-ct CROMWELL_RUNTIME] [-cm CROMWELL_MEMORY] 
                            [-cromwell CROMWELL_JAR]
```   
- `-csv`: A csv file containing information about matched tumor/normal pairs. See `MuTect2_sample.csv` for proper formatting.
- `-n`: Number of cores (default = 1).
- `-t`: Slurm job runtime. Note that this is the runtime per interval job (default = 0-9:0:0).
- `-p`: Slurm queue (default = park).
- `--mem_per_cpu`: Memory per core (default = 10G).
- `--mail_type`: Notification type (default = FAIL).
- `--mail_user`: Email.
- `-reference`: Reference FASTA file path (default = `/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta`).
- `-dbsnp`: dbSNP database VCF path (default = `/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf.gz`).
- `-data_type`: WGS or WES. 
- `-r1`: The lower range index bound for BAMs to submit from the csv file (default = 1). Index 1 is the lowest. 
- `-r2`: The upper range index bound for BAMs to submit from the csv file (default = 100000).
- `-cn`: Number of cores per Cromwell job (default = 1). 
- `-ct`: Runtime per Cromwell job in minutes (default = 1000). 
- `-cm`: Memory per Cromwell job in MB (default = 7000). 
- `-cromwell`: Jar file for Cromwell execution (default = `/n/shared_db/singularity/hmsrc-gatk/cromwell-43.jar`).

5. `HaplotypeCaller_read.py`: Wrapper script to run the GATK HaplotypeCaller pipeline for germline mutation calling. 
```
usage: python3 HaplotypeCaller_read.py [-csv TUMOR/NORMAL_CSV] [-normal NORMAL_DIRECTORY] [-output_path OUTPUT_PATH] [-p QUEUE] 
                                       [-t RUNTIME] [-r1] [-r2] [--mem_per_cpu MEM_PER_CPU] [--mail_type MAIL_TYPE] 
                                       [--mail_user MAIL_USER] [-gatk PATH_TO_GATK] [-scatter COUNT] [-reference REFERENCE.FASTA] 
                                       [-reference_name REFERENCE] [-n NUM_CORES] [-cn CROMWELL_CORES] [-ct CROMWELL_RUNTIME] 
                                       [-cm CROMWELL_MEMORY] [-cromwell CROMWELL_JAR] [-picard PICARD_PATH]  
``` 
- `-csv`: A csv file containing information about matched tumor/normal pairs. See `MuTect2_sample.csv` for proper formatting.
- `-n`: Number of cores (default = 1).
- `-t`: Slurm job runtime. Note that this is the runtime per interval job (default = 0-12:0:0).
- `-p`: Slurm queue (default = park).
- `--mem_per_cpu`: Memory per core (default = 10G).
- `--mail_type`: Notification type (default = FAIL)
- `--mail_user`: Email.
- `-gatk`: GATK executable path (default = `/home/mk446/BiO/Install/GATK4.1.2.0/gatk`).
- `-scatter`: Number of interval files (splits calling into genomic intervals to speed computation) (default = 50).
- `-r1`: The lower range index bound for BAMs to submit from the csv file (default = 0). Index 0 is the lowest. 
- `-r2`: The upper range index bound for BAMs to submit from the csv file (default = 100000).
- `-reference`: Path to reference FASTA (default = `/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta`).
- `-reference`: Reference name (b37, hg19, etc.) (default = `b37`).
- `-cn`: Number of cores per Cromwell job (default = 1). 
- `-ct`: Runtime per Cromwell job in minutes (default = 1000). 
- `-cm`: Memory per Cromwell job in MB (default = 7000). 
- `-cromwell`: Jar file for Cromwell execution (default = `/n/shared_db/singularity/hmsrc-gatk/cromwell-43.jar`).
- `-picard`: Jar file for Picard software (default = `/home/mk446/BiO/Install/picard-tools-2.5.0/picard.jar`).

6. `Intersect.sh`: Bash script to organize and intersect the calls by MuTecT and MuSE. 
```
usage: sh Intersect.sh [OUTPUT_DIRECTORY] [MUTECT2_PATH] [MUSE_PATH]
``` 
- Both the MuTecT2 and MuSE paths should be paths to the list of files directly outputted by MuTecT2 and MuSE. The script will create and organize and manipulate files on its own. 
- The MuSE path is optional, but recommended.  

7. `Filter.sh`: Bash script to filter the intersection of the calls. 
```
usage: sh Filter.sh [PATH_TO_INTERSECTION] [NORMAL] [MATCHED_NORMAL] [CSV] [ALT_CUT] [TOTAL_CUT] [VAF_CUT] [MAF_CUT]
                    [REFERENCE] [ANNOVAR_DATABASES] [BAM_PATH] [ANNOVAR_SCRIPT] [FILTER_WITH_PANEL] [PANEL]
```
- **All fields are required unless indicated. All paths should be full paths.**
- `[PATH_TO_INTERSECTION]`: The full path to the directory of the intersection of the calls. 
- `[NORMAL]`: The full path to the directory of the normal calls from HaplotypeCaller or a Panel of Normals (ie. a sample of germline calls to filter out). 
- `[MATCHED_NORMAL]`: A Boolean value indicating whether or not the normal is a matched normal (ie. from GenotypeGVCFs). 
- `[CSV]`: Path to the original csv file containing matched tumor/normal pairs. 
- `[ALT_CUT]/[VAF_CUT]`: The alternate read-level depth/VAF to cut at. These will be filtered into a file with `bad_somatic_quality` in the filename. 
- `[TOTAL_CUT]`: The total read-level depth (ie. alt + ref) to cut at. 
- `[MAF_CUT]`: The population germline cutoff to cut at. 
- `[REFERENCE]`: hg19 or hg38 (for Annovar). 
- `[ANNOVAR_DATABASES]`: The path to the Annovar databases created from `SetupDatabases.sh`. 
- `[BAM_PATH]`: The full path to the directory of BAM files.
- `[ANNOVAR_SCRIPT]`: The path to the Annovar Perl script. On Orchestra, this is `/home/mk446/bin/annovar/table_annovar.pl`. 
- `[FILTER_WITH_PANEL]`: True (if PoN filtering is desired), False (otherwise). Currently, panel filtering is only supported for hg19/b37.
- `[PANEL]` (optional): The path to a Panel of Normals to filter with, if desired. For hg19/b37, the TCGA panel located at `/n/data1/hms/dbmi/park/victor/references/` is recommended.

8. `Annotate.sh`: Bash script to annotate the filtering results and merge them into final annotated callsets. 
```
usage: sh Annotate.sh [ANNOVAR_SCRIPT] [ANNOVAR_DATABASES] [OUTPUT_DIRECTORY] [PATH_TO_MUTECT2] [REFERENCE] [CSV] [PATH_TO_NORMAL]
```
- All paths should be full paths.
- `[ANNOVAR_SCRIPT]`: The path to the Annovar Perl script. On Orchestra, this is `/home/mk446/bin/annovar/table_annovar.pl`. 
- `[ANNOVAR_DATABASES]`: The path to the Annovar databases created from `SetupDatabases.sh`. 
- `[OUTPUT_DIRECTORY]`: The same output directory used before.
- `[PATH_TO_MUTECT]`: Path to MuTecT output.
- `[REFERENCE]`: hg19 or hg38 (for Annovar). 
- `[CSV]`: Path to the original csv file containing matched tumor/normal pairs. 
- `[PATH_TO_NORMAL]` (optional): The full path to the directory of the normal calls from HaplotypeCaller (if used). 

9. `runAll.sh`: Bash script to execute all filtering steps. **NOTE: This feature is not yet available. Please follow the steps in the instructions instead**.
```
usage: sh runAll.sh [CONFIG_PATH] 
```
- `[CONFIG_PATH]`: Path to JSON configuration file.

10. `cleanFiles.py`: A file used to create the proper directory/file structure for the filtering portion of the pipeline. Instead of moving files, this will read the input csv file and create symbolic links. 
 ```
 usage: python3 cleanFiles.py [-mutect_path MUTECT_OUTPUT_PATH] [-muse_path MUSE_OUTPUT_PATH] 
                              [-haplotypecaller_path HAPLOTYPECALLER_OUTPUT_PATH] [-bam_path BAM_PATH] [-csv CSV_OF_FILES] 
                              [-out OUTPUT_DIRECTORY]
 ```
 - `-bam_path`: The path to the BAM files and their respective index files. 
 - `-csv`: The csv detailing the organization of samples. See `/postProcessing/cleanUp.csv` for proper formatting.
 - `-out`: The output directory for the tumor-normal matched csv to be written to. 

