"""
Test: 
python3 /n/data1/hms/dbmi/park/victor/scripts/other/GenomicsDBImport_read.py -in_dir /n/data1/bch/genetics/walsh-park/data/Glia_single_cells/combined_bam_list/.HaplotypeCaller -out_dir /n/data1/hms/dbmi/park/victor/Craig.2


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
    parser.add_argument('-tumor_name', '--input_tumor_name', help='name of tumor file')
    parser.add_argument('-in_dir', '--input_directory', help='input directory containing gVCFs to be merged')
    parser.add_argument('-out_dir', '--output_directory', default='./', help='directory to which the output will be written to')
    parser.add_argument('-n', '--num_cores', default='15', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='10-00:00:00', help='slurm job submission option')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='25G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='ALL', help='slurm job submission option')
    parser.add_argument('--mail_user', default='victor_mao@hms.harvard.edu', help='slurm job submission option')
    parser.add_argument('-gatk', '--gatk_path', default='/home/mk446/BiO/Install/GATK4.1.2.0/gatk', help='path to software')
    parser.add_argument('-reference', '--reference_path', default='/n/data1/hms/dbmi/park/victor/references/human_g1k_v37.fasta', help='path to reference_path file')
    parser.add_argument('-scatter', '--scatter_size', default='1')
    return parser.parse_args()

def main():

    os.system('module load gcc/6.2.0 python/3.6.0 java perl')

    args = parse_args()
    output_dir = clean_arg_paths(args)

    generate_regions_files(args)
    region_files = return_region_files(args)

    for region_file in region_files:
        slurm_command = return_slurm_command(args)
        output_file_name = gen_output_file_name_GenomicsDB(args, output_dir)

        genomics_DB_dir = os.path.join(output_dir, 'GenomicsDB')
        #normals_dir = os.path.join(args.output_directory, 'normals')
        #os.makedirs(genomics_DB_dir, exist_ok=True)
        #normal_list = collect_normals(normals_dir)
        normal_list = collect_normals(args.input_directory)

        primary_command = return_primary_command_GenomicsDB(args, region_file, genomics_DB_dir, normal_list)

        sh_file_name = gen_sh_file_name(args, output_file_name)
        write_out(args, slurm_command, primary_command, sh_file_name)

        submit_job(sh_file_name)

def clean_arg_paths(args):
    """Modifies all user-inputted directory paths such that they end with a '/'"""
    d = vars(args)
    for arg in d.keys():
        if 'output_directory' in arg and d[arg]=='./': d[arg] = os.getcwd()
        if 'directory' in arg and d[arg] is not None and d[arg][-1] != '/': d[arg] += '/'

    output_dir = os.path.join(args.output_directory, '.GenotypeGVCFs')
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

def generate_regions_files(args, output_directory):
    os.makedirs(os.path.dirname(os.path.join(output_directory, '.Mutect2/.regions/')), exist_ok=True)
    os.system(args.gatk_path + ' SplitIntervals' + '\\' + '\n' + \
     '\t' + '-R ' + args.reference_path + ' \\' + '\n' + \
     '\t' + '-scatter ' + args.scatter_size + ' \\' + '\n' + \
     '\t' + '-O ' + os.path.join(output_directory, '.Mutect2/.regions/') + ' \\')
    sleep(15)

def return_region_files(args, output_directory):
    region_files = [os.path.join(output_directory,'.Mutect2/.regions/') + file for file in os.listdir(os.path.join(output_directory, '.Mutect2/.regions/'))]
    return region_files

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
    output_file_name = os.path.join(output_directory, 'dbCreate/GenomicsDB')
    return output_file_name

def return_primary_command_GenomicsDB(args, region_file, db_path, normal_path_list):
    primary_command = args.gatk_path + ' GenomicsDBImport' + \
    " --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true'" + \
    ' -R ' + args.reference_path + \
    ' -L ' + region_file + \
    ' --genomicsdb-workspace-path ' + db_path
    for normal in normal_path_list:
        primary_command += ' -V ' + normal 

    return primary_command

def collect_normals(normals_path):
    normals = [os.path.realpath(file) for file in glob.glob(os.path.join(normals_path, '*.gz'))]
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

if __name__ == "__main__":
    main()



    
