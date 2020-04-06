# SNVCurate
A filtering pipeline for producing a callset of somatic mutations. It also provides a method to curate targeted capture regions in the event that this BED file is not available. 
This pipeline is used only for running on Orchestra, the Harvard Medical School cluster. 

**Note: This pipeline presumes that the BAM files are formatted properly (demultiplexed and indexed correctly). If they are not, please see BAMCurate first before running this pipeline. If no matched normal is available, run this pipeline using another normal of the same sequencing type.**

## Running the SNV curating pipeline: 
### Environment/file setup:
1. Create two files: One csv with tumor/normal matched pairs (leave `N/A` if no normal) and one text file with each normal sample (if applicable). See `Tumor_Normal_sample.csv` and `HaplotypeCaller_sample.txt` for formatting.
2. Set up the environment using Conda, or download and install the packages individually (see below for dependencies). 

### Running the pipeline (together): 
There are two ways to run this pipeline: 
1. Create a config file `parameters.config`. See `parameters_example.config` for details.
2. Run `runAll.sh`.

### Running the pipeline (individual): 
1. Run MuTecT2 on the BAM files using `Mutect2_read.py`. 
2. If you desire higher specifity for somatic calls (recommended), run MuSE using `MuSE_read.py`. 
3. If there are matched normals, run HaplotypeCaller on these normals using `HaplotypeCaller_read.py`.
4. Run the script `Clean.sh` to clean the calls made by the callers. 
5. Run `GenotypeGVCFs_read.py` on the new HaplotypeCaller gVCFs to finish calling germline mutations.
4. Run the script `Intersect.sh` to intersect the two calls. 
5. Run the script `Filter.sh` to filter the intersection.
6. Run the script `Annotate.sh` to finish. 

## Dependencies: 
- _libgcc_mutex=0.1=conda_forge
- _openmp_mutex=4.5=0_gnu
- bzip2=1.0.8=h516909a_2
- ca-certificates=2019.11.28=hecc5488_0
- certifi=2019.11.28=py37hc8dfbb8_1
- curl=7.68.0=hf8cf82a_0
- jq=1.6=h14c3975_1000
- krb5=1.16.4=h2fd8d38_0
- ld_impl_linux-64=2.34=h53a641e_0
- libblas=3.8.0=14_openblas
- libcblas=3.8.0=14_openblas
- libcurl=7.68.0=hda55be3_0
- libdeflate=1.5=h516909a_0
- libedit=3.1.20170329=hf8c457e_1001
- libffi=3.2.1=he1b5a44_1007
- libgcc-ng=9.2.0=h24d8f2e_2
- libgfortran-ng=7.3.0=hdf63c60_5
- libgomp=9.2.0=h24d8f2e_2
- liblapack=3.8.0=14_openblas
- libopenblas=0.3.7=h5ec1e0e_6
- libssh2=1.8.2=h22169c7_2
- libstdcxx-ng=9.2.0=hdf63c60_2
- muse=1.0.rc=hdbcaa40_3
- ncurses=6.1=hf484d3e_1002
- numpy=1.18.1=py37h8960a57_1
- oniguruma=6.9.3=h516909a_0
- openssl=1.1.1f=h516909a_0
- pandas=1.0.3=py37h0da4684_0
- pip=20.0.2=py_2
- pysam=0.15.4=py37hbcae180_0
- python=3.7.6=h8356626_5_cpython
- python-dateutil=2.8.1=py_0
- python_abi=3.7=1_cp37m
- pytz=2019.3=py_0
- readline=8.0=hf8c457e_0
- setuptools=46.1.3=py37hc8dfbb8_0
- six=1.14.0=py_1
- sqlite=3.30.1=hcee41ef_0
- tabix=0.2.6=ha92aebf_0
- tk=8.6.10=hed695b0_0
- wheel=0.34.2=py_1
- xz=5.2.4=h516909a_1002
- zlib=1.2.11=h516909a_1006

## Information about relevant scripts: 
1. `Mutect2_read.py`: Wrapper script to run the GATK MuTect2 pipeline for somatic mutation calling. 
```
usage: python3 Mutect2_read.py [-tumor INPUT_DIRECTORY] [-normal NORMAL_DIRECTORY] [-out OUTPUT_DIRECTORY] [-csv TUMOR/NORMAL_CSV] 
                     [-pon PANEL_OF_NORMALS] [-n NUM_CORES] [-t RUNTIME] [-p QUEUE] [--mem_per_cpu MEM_PER_CPU] [--mail_type MAIL_TYPE]
                     [--mail_user MAIL_USER] [-gatk4 PATH_TO_GATK4] [-reference REFERENCE.FASTA] [-dbsnp dbSNP.vcf] [-gnomad GNOMAD.vcf] 
                     [-scatter COUNT] [-r1 START] [-r2 END]
```
- `-pon`: Used when there are no matched normals. 
- `-csv`: A csv file containing information about matched tumor/normal pairs. See `MuTect2_sample.csv` for proper formatting.
- `-n`: Number of cores (default = 1).
- `-t`: Slurm job runtime. Note that this is the runtime per interval job (default = 0-12:0:0).
- `-p`: Slurm queue (default = park).
- `--mem_per_cpu`: Memory per core (default = 10G).
- `--mail_type`: Notification type (default = FAIL)
- `--mail_user`: Email.
- `-gatk4`: GATK executable path(default = `/n/data1/hms/dbmi/park/alon/software/gatk/gatk-4.0.3.0/gatk`).
- `-reference`: Reference FASTA file path (default = `/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta`).
- `-dbsnp`: dbSNP database VCF path (default = `/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf`).
- `-scatter`: Number of interval files (splits calling into genomic intervals to speed computation) (default = 50).
- `-r1`: The lower range index bound for BAMs to submit from the csv file (default = 1). Index 1 is the lowest. 
- `-r2`: The upper range index bound for BAMs to submit from the csv file (default = 100000).

2. `MuSE_read.py`: Wrapper script to run the MuSE pipeline for somatic mutation calling. 
```
usage: python3 MuSE_read.py [-tumor INPUT_DIRECTORY] [-normal NORMAL_DIRECTORY] [-out OUTPUT_DIRECTORY] [-csv TUMOR/NORMAL_CSV] [-n NUM_CORES]
                     [-t RUNTIME] [-p QUEUE] [--mem_per_cpu MEM_PER_CPU] [--mail_type MAIL_TYPE] [--mail_user MAIL_USER] 
                     [-reference REFERENCE.FASTA] [-dbsnp dbSNP.vcf] [-mode MODE] [-data_type WGS/WES] [-r1 START] [-r2 END]
```   
- `-csv`: A csv file containing information about matched tumor/normal pairs. See `MuTect2_sample.csv` for proper formatting.
- `-n`: Number of cores (default = 1).
- `-t`: Slurm job runtime. Note that this is the runtime per interval job (default = 0-9:0:0).
- `-p`: Slurm queue (default = park).
- `--mem_per_cpu`: Memory per core (default = 10G).
- `--mail_type`: Notification type (default = FAIL).
- `--mail_user`: Email.
- `-reference`: Reference FASTA file path (default = `/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta`).
- `-dbsnp`: dbSNP database VCF path (default = `/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf`).
- `-mode`: Either `call` or `sump`. You must call before running `sump`. Note that in `sump`, you do not need a normal file and instead should input the resulting file from `call` using the `-tumor` flag. 
- `-data_type`: WGS or WES. 
- `-r1`: The lower range index bound for BAMs to submit from the csv file (default = 1). Index 1 is the lowest. 
- `-r2`: The upper range index bound for BAMs to submit from the csv file (default = 100000).

3. `HaplotypeCaller_read.py`: Wrapper script to run the GATK HaplotypeCaller pipeline for germline mutation calling. 
```
usage: python3 HaplotypeCaller_read.py [-input_path INPUT_PATH] [-output_path OUTPUT_PATH] [-p QUEUE] [-t RUNTIME] [-r1] [-r2] 
                                  [--mem_per_cpu MEM_PER_CPU] [--mail_type MAIL_TYPE] [--mail_user MAIL_USER] [-gatk PATH_TO_GATK]
                                    [-scatter COUNT] [-reference REFERENCE.FASTA] [-n NUM_CORES]
``` 
- `-input_path`: Path to a text file that lists the full path to each matched normal. See `HaplotypeCaller_sample.txt` as an example.
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

4. `GenotypeGVCFs_read.py`: Wrapper script to run the GATK GenotypeGVCFs pipeline for germline mutation calling. 
```
usage: python3 GenotypeGVCFs_read.py [-in_dir INPUT_PATH] [-out_dir OUTPUT_PATH] [-p QUEUE] [-t RUNTIME]                                                                        [--mem_per_cpu MEM_PER_CPU] [--mail_type MAIL_TYPE] [--mail_user MAIL_USER] [-gatk PATH_TO_GATK]
                                       [-reference REFERENCE.FASTA] [-n NUM_CORES]
```
- `-in_dir`: Path to a text file that lists the full path to each matched normal. See `HaplotypeCaller_sample.txt` as an example.
- `-n`: Number of cores (default = 1).
- `-t`: Slurm job runtime. Note that this is the runtime per interval job (default = 0-12:0:0).
- `-p`: Slurm queue (default = park).
- `--mem_per_cpu`: Memory per core (default = 10G).
- `--mail_type`: Notification type (default = FAIL)
- `--mail_user`: Email.
- `-gatk`: GATK executable path (default = `/home/mk446/BiO/Install/GATK4.1.2.0/gatk`).
- `-reference`: Path to reference FASTA (default = `/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta`). 
- `-mode`: Genotyping in joint/group mode or single sample mode (default = `single`). 


5. `Clean.sh`: Bash script to organize, merge, and clean files outputted by GATK callers. **Note: This must be run after MuTecT, MuSE, and HaplotypeCaller but before GenotypeGVCFs.** 
```
usage: sh Clean.sh [GENOTYPEGVCFs_PATH] [MUTECT2_PATH] [MUSE_PATH]
``` 
- If there is no matched normal, then input `none` for `[GENOTYPEGVCFs_PATH]`.
- The MuSE path is optional, but recommended. 


6. `Intersect.sh`: Bash script to organize and intersect the calls by MuTecT and MuSE. 
```
usage: sh Intersect.sh [OUTPUT_DIRECTORY] [MUTECT2_PATH] [MULTI-LANE] [MUSE_PATH]
``` 
- Both the MuTecT2 and MuSE paths should be paths to the list of files directly outputted by MuTecT2 and MuSE. The script will create and organize and manipulate files on its own. 
- The MuSE path is optional, but recommended. 
- `[MULTI-LANE]` is a Boolean referring to if the somatic callset by MuTecT has multiple lanes/sample IDs. The script will merge them into one before filtering, thus still preserving all calls. 

7. `Filter.sh`: Bash script to filter the intersection of the calls. 
```
usage: sh Filter.sh [PATH_TO_INTERSECTION] [NORMAL] [MATCHED_NORMAL] [CSV] [PANEL] [ALT_CUT] [TOTAL_CUT] [VAF_CUT] [MAF_CUT]                                   [REFERENCE] [PATH_TO_ANNOVAR_DATABASES] [BAM_PATH] [FILTER_WITH_PANEL] [PANEL]
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
- `[PATH_TO_ANNOVAR_DATABASES]`: A path to a directory for which the script will output soft links to Annovar databases, along with custom filtering BED files. 
- `[BAM_PATH]`: The full path to the directory of BAM files.
- `[FILTER_WITH_PANEL]`: True (if PoN filtering is desired), False (otherwise). 
- `[PANEL]` (optional): The path to a Panel of Normals to filter with, if desired.

8. `Annotate.sh`: Bash script to annotate the filtering results and merge them into final annotated callsets. 
```
usage: sh Annotate.sh [OUTPUT_DIRECTORY] [PATH_TO_MUTECT2] [REFERENCE] [CSV] [PATH_TO_NORMAL]
```
Additional Information/Default parameters:  
- All paths should be full paths.
- `[OUTPUT_DIRECTORY]`: The same output directory used before.
- `[PATH_TO_MUTECT]`: Path to MuTecT output.
- `[REFERENCE]`: hg19 or hg38 (for Annovar). 
- `[CSV]`: Path to the original csv file containing matched tumor/normal pairs. 
- `[PATH_TO_NORMAL]`: The full path to the directory of the normal calls from HaplotypeCaller (if used). 

    
