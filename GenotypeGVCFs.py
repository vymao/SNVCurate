"""
Test: 
python3 /n/data1/hms/dbmi/park/victor/scripts/other/GenotypeGVCFs_read.2.py -in_dir /n/data1/bch/genetics/walsh-park/data/Glia_single_cells/combined_bam_list/.HaplotypeCaller -out_dir /n/data1/hms/dbmi/park/victor/Craig/11.19.2019_with_Alice/ -t 3-23:00:00 --mem_per_cpu 50G -n 8 -p park -reference /home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta -temp_dir /n/data1/hms/dbmi/park/victor/Craig/11.19.2019_with_Alice/

python3 /n/data1/hms/dbmi/park/victor/scripts/other/GenotypeGVCFs_read.2.py -in_dir /n/data1/hms/dbmi/park/victor/other/tests/Craigtest -out_dir /n/data1/hms/dbmi/park/victor/other/tests/Craigtest

python3 /n/data1/hms/dbmi/park/victor/scripts/other/GenotypeGVCFs_read.2.py -in_dir /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK/HaplotypeCaller_runs -out_dir /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK/GenotypeGVCFs -mode single -p medium 


python3 /n/data1/hms/dbmi/park/victor/scripts/other/GenotypeGVCFs_read.2.py -in_dir /n/data1/hms/dbmi/park/victor/Doga/Gerburg_WES_2/normal_list.HaplotypeCaller/.HaplotypeCaller -out_dir /n/data1/hms/dbmi/park/victor/Doga/Gerburg_WES_2/normal_list.HaplotypeCaller/GenotypeGVCFs -mode single -p park


python3 /n/data1/hms/dbmi/park/victor/scripts/other/GenotypeGVCFs_read.2.py -in_dir /home/clb36/parkhome/Juan/original/bam_list.HaplotypeCaller/HaplotypeCaller_merged -out_dir /home/clb36/parkhome/Juan/original/bam_list.HaplotypeCaller/GenotypeGVCFs -t 5-23:00:00 --mem_per_cpu 50G -n 8 -p park -reference /n/data1/hms/dbmi/park/victor/references/mouse/ncbi-genomes-2019-11-21/GCF_000001635.26_GRCm38.p6_genomic.fasta -temp_dir /home/clb36/parkhome/Juan/original/bam_list.HaplotypeCaller/GenotypeGVCFs


python3 /n/data1/hms/dbmi/park/victor/scripts/other/GenotypeGVCFs_read.2.py -in_dir /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/HaplotypeCaller_recalibrated/normal_list.HaplotypeCaller/.HaplotypeCaller/new -out_dir /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/HaplotypeCaller_recalibrated/GenotypeGVCFs -mode single -p park

"""



import argparse
import os
import re
import ntpath
import glob
import sys
from time import sleep

def parse_args():
    """Uses argparse to enable user to customize script functionality"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-in_dir', '--input_directory', help='input directory containing gVCFs to be merged')
    parser.add_argument('-out_dir', '--output_directory', default='./', help='directory to which the output will be written to')
    parser.add_argument('-n', '--num_cores', default='2', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='2-00:00:00', help='slurm job submission option')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='50G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', default='victor_mao@hms.harvard.edu', help='slurm job submission option')
    parser.add_argument('-gatk3', '--gatk3_path', default='/home/mk446/BiO/Install/GATK3.8-0/GenomeAnalysisTK.jar', help='path to GATK3 jar script')
    parser.add_argument('-reference', '--reference_path', default='/n/data1/hms/dbmi/park/victor/references/human_g1k_v37.fasta', help='path to reference_path file')
    parser.add_argument('-mode', default='group', help='Call Joint Genotyping in group mode or single sample mode')
    parser.add_argument('-r1', default=0, help='Lower index bound for single sample mode')
    parser.add_argument('-r2', default=1000000, help='Upper index bound for single sample mode')
    parser.add_argument('-temp_dir', default=None, help='Temporary directory for Java')
    return parser.parse_args()

def clean_arg_paths(args):
    """Modifies all user-inputted directory paths such that they end with a '/'"""
    d = vars(args)
    for arg in d.keys():
        if 'output_directory' in arg and d[arg]=='./': d[arg] = os.getcwd()
        if 'directory' in arg and d[arg] is not None and d[arg][-1] != '/': d[arg] += '/'

def return_slurm_command(args):
    """Returns slurm command given args provided"""
    slurm_command = '#!/bin/bash\n' + \
                '#SBATCH -n ' + args.num_cores + '\n' + \
                '#SBATCH -t ' + args.runtime + '\n' + \
                '#SBATCH -p ' + args.queue + '\n' + \
                '#SBATCH --mem-per-cpu=' + args.mem_per_cpu + '\n' + \
                '#SBATCH --mail-type=' + args.mail_type + '\n' + \
                '#SBATCH --mail-user=' + args.mail_user + '\n' + \
                '#SBATCH -o ' + os.path.join(args.output_directory, 'slurm_log') + '\n' + \
                '#SBATCH -e ' + os.path.join(args.output_directory, 'slurm_err') + '\n'
    if args.queue in ['park', 'priopark']:
        slurm_command += '#SBATCH --account=park_contrib' + '\n'
    slurm_command += 'module load java/jdk-1.8u112' + '\n'
    slurm_command += 'ulimit -c unlimited' + '\n'
    slurm_command += 'ulimit -s unlimited' + '\n'
    return slurm_command

def gen_output_file_name_GenomicsDB(args):
    output_file_name = os.path.join(args.output_directory, 'GenotypeGVCFs')
    return output_file_name

def gen_sample_output_file_name_GenomicsDB(args, sample_name):
    sample = ntpath.basename(sample_name).rsplit('.g.vcf.gz')[0]
    output_file_name = os.path.join(args.output_directory, sample + '.vcf')
    return output_file_name

def return_primary_command_GVCF(args, normal_path_list, output_file_name):
    primary_command = 'java -Xmx200g' 
    if args.temp_dir is not None: 
        primary_command += ' -Djava.io.tmpdir=' + args.temp_dir
    primary_command +=  ' -jar ' + args.gatk3_path + ' -T GenotypeGVCFs' + \
    ' -R ' + args.reference_path + \
    ' -o ' + output_file_name + \
    ' -nt ' + args.num_cores
    for normal in normal_path_list:
        primary_command += ' --variant ' + normal 

    return primary_command

def collect_normals(normals_path):
    normals = [os.path.realpath(file) for file in glob.glob(os.path.join(normals_path, '*.g.vcf.gz'))]
    return normals

def gen_sh_file_name(args, output_file_name):
    """Generates sh file name"""
    sh_file_name = os.path.dirname(output_file_name) + '/.sh/' + os.path.basename(output_file_name) + '.sh'
    return sh_file_name

def write_out(args, slurm_command, primary_command, sh_file_name):
    """"""
    os.makedirs(os.path.dirname(sh_file_name), exist_ok=True)
    with open(sh_file_name, 'w') as file:
        file.write(slurm_command + primary_command)

def submit_job(sh_file_name):
    os.chdir(os.path.dirname(sh_file_name))
    os.system('chmod +x ' + os.path.basename(sh_file_name))
    os.system('sbatch ./' + os.path.basename(sh_file_name))

def main():

    os.system('module load gcc/6.2.0 python/3.6.0 java perl')

    args = parse_args()
    clean_arg_paths(args)

    if args.mode.lower() == 'group':
        slurm_command = return_slurm_command(args)
        output_file_name = gen_output_file_name_GenomicsDB(args)

        normal_list = collect_normals(args.input_directory)

        primary_command = return_primary_command_GVCF(args, normal_list, output_file_name)

        sh_file_name = gen_sh_file_name(args, output_file_name)
        write_out(args, slurm_command, primary_command, sh_file_name)

        submit_job(sh_file_name)
        #print(primary_command)
    elif args.mode.lower() == 'single':
        slurm_command = return_slurm_command(args)
        normal_list = collect_normals(args.input_directory)

        for file in normal_list:
            if normal_list.index(file) in range(int(args.r1), int(args.r2)):
                output_file_name = gen_sample_output_file_name_GenomicsDB(args, file)
                primary_command = return_primary_command_GVCF(args, [file], output_file_name)

                sh_file_name = gen_sh_file_name(args, output_file_name)
                write_out(args, slurm_command, primary_command, sh_file_name)

                if not os.path.exists(output_file_name):
                    submit_job(sh_file_name)

if __name__ == "__main__":
    main()



    
