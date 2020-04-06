"""
Test: 
python3 /n/data1/hms/dbmi/park/victor/scripts/other/GenotypeGVCFS_read.py -in_dir /n/data1/bch/genetics/walsh-park/data/Glia_single_cells/combined_bam_list/.HaplotypeCaller -out_dir /n/data1/hms/dbmi/park/victor/Craig


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
    parser.add_argument('-in_dir', '--input_directory', help='input directory containing GenomicsDB or HaplotypeCaller output')
    parser.add_argument('-out_dir', '--output_directory', default='./', help='directory to which the output will be written to')
    parser.add_argument('-n', '--num_cores', default='1', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='0-12:00:00', help='slurm job submission option')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='10G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', help='slurm job submission option')
    parser.add_argument('-gatk', '--gatk_path', default='/home/mk446/BiO/Install/GATK4.1.2.0/gatk', help='path to GATK4 software')
    parser.add_argument('-reference', '--reference_path', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta', help='path to reference_path file')
    parser.add_argument('-mode', default='single', help='Call Joint Genotyping in group mode or single sample mode')
    return parser.parse_args()

def main():
    os.system('module load gcc/6.2.0 python/3.6.0 java perl')

    args = parse_args()
    output_dir = clean_arg_paths(args)
    if args.mode.lower() == 'group':
        slurm_command = return_slurm_command(args)
        output_file_name = gen_output_file_name_GenomicsDB(args, output_dir)
        primary_command = return_primary_command_GenomicsDB(args, output_file_name)

        sh_file_name = gen_sh_file_name(args, output_file_name)
        write_out(args, slurm_command, primary_command, sh_file_name)
     
        submit_job(sh_file_name)
    elif args.mode.lower() == 'single': 
        slurm_command = return_slurm_command(args)
        normal_list = collect_normals(args.input_directory)
        #print(normal_list)
        for file in normal_list:
            sample = os.path.basename(file)
            output_file_name = gen_output_file_name_single(file, output_dir)
            primary_command = return_primary_command_single(args, file, output_file_name)

            sh_file_name = gen_sh_file_name(args, output_file_name)
            write_out(args, slurm_command, primary_command, sh_file_name)

            if not os.path.exists(output_file_name):
                submit_job(sh_file_name)
    else: 
        print('Invalid mode.')

def clean_arg_paths(args):
    """Modifies all user-inputted directory paths such that they end with a '/'"""
    d = vars(args)
    for arg in d.keys():
        if 'output_directory' in arg and d[arg]=='./': d[arg] = os.getcwd()
        if 'directory' in arg and d[arg] is not None and d[arg][-1] != '/': d[arg] += '/'

    output_dir = os.path.join(args.output_directory, '.GenotypeGVCFs')
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def return_slurm_command(args):
    """Returns slurm command given args provided"""
    slurm_command = '#!/bin/bash\n' + \
                '#SBATCH -n ' + args.num_cores + '\n' + \
                '#SBATCH -t ' + args.runtime + '\n' + \
                '#SBATCH -p ' + args.queue + '\n' + \
                '#SBATCH --mem-per-cpu=' + args.mem_per_cpu + '\n' + \
                '#SBATCH --mail-type=' + args.mail_type + '\n' + \
                '#SBATCH --mail-user=' + args.mail_user + '\n'
    if args.queue in ['park', 'priopark']:
        slurm_command += '#SBATCH --account=park_contrib' + '\n'
    slurm_command += 'module load java/jdk-1.8u112' + '\n'
    return slurm_command

def gen_output_file_name_GenomicsDB(args, output_directory):
    output_file_name = os.path.join(output_directory, 'Genotype_grouped.vcf')
    return output_file_name

def gen_output_file_name_single(sample, output_directory):
    base = os.path.basename(sample).split('.')[0]    
    output_file_name = os.path.join(output_directory, base + '.vcf')
    return output_file_name

def return_primary_command_GenomicsDB(args, output):
    primary_command = args.gatk_path + \
    ' --java-options "-Xmx200g" GenotypeGVCFs' + \
    ' -R ' + args.reference_path + \
    ' -V ' + 'gendb://' + args.input_directory + \
    ' -O ' + output
    return primary_command

def return_primary_command_single(args, sample, output):
    primary_command = args.gatk_path + \
    ' --java-options "-Xmx200g" GenotypeGVCFs' + \
    ' -R ' + args.reference_path + \
    ' -V ' + sample + \
    ' -O ' + output 
    return primary_command

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

def collect_normals(normals_path):
    normals = [os.path.realpath(file) for file in glob.glob(os.path.join(normals_path, '*.g.vcf.gz'))]
    return normals

if __name__ == "__main__":
    main()



    
