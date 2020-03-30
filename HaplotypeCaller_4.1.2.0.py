"""
python3 /n/data1/hms/dbmi/park/victor/scripts/GATK_Germline_SNPs_Indels/HaplotypeCaller_4.1.2.0.py -tumor /n/data1/hms/dbmi/park/ethan/GERBURG/.PreProcessing/Sample_DS-bkm-085-N.txt.b37.bam -out /n/data1/hms/dbmi/park/victor/Doga/test/Gerburg_WES_test 



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
    parser.add_argument('-in_file', '--input_tumor_path', help='path to input tumor file')
    parser.add_argument('-out', '--output_directory', default='./', help='directory to which the output directory "/.HaplotypeCaller/" will be written to')
    parser.add_argument('-n', '--num_cores', default='1', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='2-12:00:00', help='slurm job submission option')
    parser.add_argument('-t2', help='dummy')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='20G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', help='slurm job submission option')
    parser.add_argument('-gatk', '--gatk_path', default='/home/mk446/BiO/Install/GATK4.1.2.0/gatk', help='path to software execution script')
    #parser.add_argument('-reference', '--reference_path', default='/n/data1/hms/dbmi/park/victor/references/human_g1k_v37.fasta', help='path to reference_path file')
    #parser.add_argument('-reference', '--reference_path', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta', help='path to reference_path file')
    parser.add_argument('-reference', '--reference_path', default='/n/data1/hms/dbmi/park/victor/references/mouse/ncbi-genomes-2019-11-21/GCF_000001635.26_GRCm38.p6_genomic.fasta', help='path to reference_path file')
    parser.add_argument('-scatter', '--scatter_size', default='50')
    return parser.parse_args()

def clean_arg_paths(args):
    """Modifies all user-inputted directory paths such that they end with a '/'"""
    d = vars(args)
    for arg in d.keys():
        if 'output_directory' in arg and d[arg]=='./': d[arg] = os.getcwd() + '/.' + ntpath.basename(args.input_tumor_path) + '_v_' + ntpath.basename(args.input_normal_path)
        if 'directory' in arg and d[arg] is not None and d[arg][-1] != '/': d[arg] += '/'

def generate_regions_files(args):
    os.makedirs(os.path.dirname(args.output_directory + '.HaplotypeCaller/.regions/'), exist_ok=True)
    os.system(args.gatk_path + ' SplitIntervals' + '\\' + '\n' + \
     '\t' + '-R ' + args.reference_path + ' \\' + '\n' + \
     '\t' + '-scatter ' + args.scatter_size + ' \\' + '\n' + \
     '\t' + '-O ' + args.output_directory + '.HaplotypeCaller/.regions/' + ' \\')
    sleep(15)

def return_region_files(args):
    region_files = [args.output_directory + '.HaplotypeCaller/.regions/' + file for file in os.listdir(args.output_directory + '.HaplotypeCaller/.regions/')]
    return region_files

def return_slurm_command(args):
    """Returns slurm command given args provided"""

    basename = ntpath.basename(args.input_tumor_path)
    log_file = '.' + basename + '.Log'
    error_file = '.' + basename + '.Errors'
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

def gen_output_file_name(args, region_file):
    output_file_name = args.output_directory + '.HaplotypeCaller/' + ntpath.basename(args.input_tumor_path) + '_' + ntpath.basename(region_file[1:]) + '.g.vcf'

    return output_file_name

def return_primary_command(args, output_file_name, region_file):
    primary_command = args.gatk_path + ' --java-options "-Xmx100g" HaplotypeCaller ' + \
     '-R ' + args.reference_path + \
     ' -I ' + args.input_tumor_path  + \
     ' -L ' + args.output_directory + '.HaplotypeCaller/.regions/' + ntpath.basename(region_file) + \
     ' -O ' + output_file_name + \
     ' -ERC GVCF'
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

def main():
    os.system('module load gcc/6.2.0 python/3.6.0 samtools/1.3.1 bwa/0.7.15 java')  

    args = parse_args()
    clean_arg_paths(args)

    generate_regions_files(args)
    region_files = return_region_files(args)

    for region_file in region_files:
        slurm_command = return_slurm_command(args)
        output_file_name = gen_output_file_name(args, region_file)

        primary_command = return_primary_command(args, output_file_name, region_file)

        sh_file_name = gen_sh_file_name(args, output_file_name)
        write_out(args, slurm_command, primary_command, sh_file_name)

        sample_name = ntpath.basename(sh_file_name).split('.bam')[0]
        mutect_path = ntpath.dirname(ntpath.dirname(sh_file_name))
        path_to_vcf = os.path.join(mutect_path, sample_name)
        vcf = ntpath.basename(sh_file_name).split('.sh')[0]
        path_to_vcf = os.path.join(path_to_vcf, vcf)    
        #if not os.path.isfile(path_to_vcf):
        #    submit_job(sh_file_name)
            #print(primary_command)
        submit_job(sh_file_name)
if __name__ == "__main__":
    main()
