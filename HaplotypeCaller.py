
"""
Execution test:
python3 /n/data1/hms/dbmi/park/victor/scripts/GATK_Germline_SNPs_Indels/HaplotypeCaller.py -in_file /n/data1/hms/dbmi/park/victor/other/tests/5817_liver_bulk.bam -out /n/data1/hms/dbmi/park/victor/other/tests/ -t2 2100

python3 /n/data1/hms/dbmi/park/victor/scripts/GATK_Germline_SNPs_Indels/HaplotypeCaller.py -in_dir /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI 

python3 /n/data1/hms/dbmi/park/victor/scripts/GATK_Germline_SNPs_Indels/HaplotypeCaller.py -in_dir /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI 

python3 /n/data1/hms/dbmi/park/victor/scripts/GATK_Germline_SNPs_Indels/HaplotypeCaller.py -in_file s_DS_bkm_085_N_bc0069_Proj_5065_E_L000_mrg_cl_aln_srt_MD_IR_BR.bam -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI -p medium

python3 /n/data1/hms/dbmi/park/victor/scripts/GATK_Germline_SNPs_Indels/HaplotypeCaller.py -in_file /n/data1/hms/dbmi/park/victor/Doga/test/s_DS_bkm_085_N_bc0069_Proj_5065_E_L000_mrg_cl_aln_srt_MD_IR_BR.bam -out /n/data1/hms/dbmi/park/victor/Doga/test -p park



"""
from time import sleep
import argparse
import os
import re
import ntpath
import glob
import sys
import fileinput
from shutil import copyfile

def parse_args():
    """Uses argparse to enable user to customize script functionality"""
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in_dir', '--input_directory', help='path to directory containing input files')
    parser.add_argument('-in_file', '--input_file_path', help='path to input file')
    parser.add_argument('-out', '--output_directory', default='./', help='directory to which the output directory "/.HaplotypeCaller/" will be written to')
    parser.add_argument('-L', '--scattered_calling_intervals_list', default='/n/data1/hms/dbmi/park/alon/software/gatk-test-data/intervals/b37_wgs_scattered_calling_intervals.txt', help='Reference genome-dependent. Scattered calling interval list')
    parser.add_argument('-n', '--num_cores', default='4', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='10-23:59:00', help='slurm job submission option')
    parser.add_argument('-t2', '--cromruntime', default='4500', help='cromwell run time; please specify as number of minutes')
    parser.add_argument('-p', '--queue', default='long', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='15G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', default='victor_mao@hms.harvard.edu', help='slurm job submission option')
    #parser.add_argument('-overrides', '--overrides_path', default='/n/data1/hms/dbmi/park/alon/software/gatk4-data-processing-master/overrides.conf.non_park', help='path to overrides.conf file')
    #parser.add_argument('-cromwell', '--cromwell_path', default='/n/data1/hms/dbmi/park/alon/software/cromwell-36.jar', help='path to cromwell.jar file')
    parser.add_argument('-overrides', '--overrides_path', default='/n/data1/hms/dbmi/park/victor/scripts/GATK_Germline_SNPs_Indels/overrides.conf.new.non_park', help='path to overrides.conf file')
    parser.add_argument('-cromwell', '--cromwell_path', default='/n/shared_db/singularity/hmsrc-gatk/cromwell-43.jar', help='path to cromwell.jar file')

    parser.add_argument('-gatk_wdl', '--gatk4_haplotypecaller_path', default='/n/data1/hms/dbmi/park/alon/software/gatk4-germline-snps-indels-master/haplotypecaller-gvcf-gatk4_shorttime.wdl', help='path to gatk4-haplotypecaller file')
    parser.add_argument('-input_json', '--input_json_path', default='/n/data1/hms/dbmi/park/alon/software/gatk4-germline-snps-indels-master/haplotypecaller-gvcf-gatk4.hg37.wgs.inputs.json', help='Reference genome-dependent. Path to gatk4-haplotypecaller file')
    return parser.parse_args()

def clean_arg_paths(args):
    """Modifies all user-inputted directory paths such that they end with a '/'"""
    d = vars(args)
    for arg in d.keys():
        if 'input_directory' in arg and d[arg]=='./': d[arg] = os.getcwd()
        if 'output_directory' in arg and d[arg]=='./': d[arg] = os.getcwd()    
        if 'directory' in arg and d[arg] is not None and d[arg][-1] != '/': d[arg] += '/'

def return_input_files(args, ext):
    input_bams = [os.path.realpath(file) for file in glob.glob(args.input_directory + '*.' + ext)]
    return input_bams

def sort_by_size(input_files):
    for i in range(len(input_files)):
        input_files[i] = (input_files[i], os.path.getsize(input_files[i]))
    input_files.sort(key=lambda filename: filename[1])
    for i in range(len(input_files)):
        input_files[i] = input_files[i][0]
    return input_files

def generate_input_json(args, input_file):
    dir = args.output_directory + '.HaplotypeCaller/' + '.' + os.path.basename(input_file) + '/'
    os.makedirs(dir, exist_ok=True)

    bai_suffix = '.bai'
    path = re.sub('.bam', '.bam.bai', input_file)
    if os.path.isfile(path) and os.access(path, os.R_OK):
        bai_suffix = '.bam.bai'
    
    copyfile(args.input_json_path, dir + 'input.json')

    with open(dir + 'input.json', "r") as f:
        for _ in range(13):
            line = f.readline()
    string = line.split(':')[1].split(',')[0].strip(' ').strip('"')

    with fileinput.FileInput(dir + 'input.json', inplace=True) as file:
        for line in file:
            print(line.replace(
            "gs://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bam", input_file).replace(
            "gs://gatk-test-data/wgs_bam/NA12878_24RG_hg38/NA12878_24RG_small.hg38.bai", re.sub('.bam', bai_suffix, input_file)).replace(
            "output_dir/", args.output_directory + '.HaplotypeCaller/').replace(
            string, args.scattered_calling_intervals_list), end='')

    copyfile(args.overrides_path, dir + 'overrides.conf')
    with fileinput.FileInput(dir + 'overrides.conf', inplace=True) as file:
        for line in file:
            print(line.replace(
            "medium", args.queue), end='')
    with fileinput.FileInput(dir + 'overrides.conf', inplace=True) as file:
        for line in file:
            print(re.sub(r"runtime_minutes = [0-9]*", 'runtime_minutes = ' + args.cromruntime, line), end='')

    copyfile(args.gatk4_haplotypecaller_path, dir + 'haplotypecaller-gvcf-gatk4.wdl')
    with fileinput.FileInput(dir + 'haplotypecaller-gvcf-gatk4.wdl', inplace=True) as file:
        for line in file:
            print(line.replace(
            "priopark", args.queue), end='')
    with fileinput.FileInput(dir + 'haplotypecaller-gvcf-gatk4.wdl', inplace=True) as file:
        for line in file: 
            print(re.sub(r'runtime_minutes: \".*\"', 'runtime_minutes: ' + args.cromruntime, line), end='')

    return dir + 'input.json'

def return_slurm_command(args):
    """Returns slurm command given args provided"""
    slurm_command = '#!/bin/bash\n' + \
                '#SBATCH -n ' + args.num_cores + '\n' + \
                '#SBATCH -t ' + args.runtime + '\n' + \
                '#SBATCH -p ' + args.queue + '\n' + \
                '#SBATCH --mem-per-cpu=' + args.mem_per_cpu + '\n' + \
                '#SBATCH --mail-type=' + args.mail_type + '\n' + \
                '#SBATCH --mail-user=' + args.mail_user + '\n' + \
        '#SBATCH --exclude=compute-p-17-[34-46]' + '\n'
    if args.queue in ['park', 'priopark']:
        slurm_command += '#SBATCH --account=park_contrib' + '\n'
    return slurm_command

def gen_output_file_name(args, input_file):
    output_file_name = args.output_directory + '.HaplotypeCaller/' + '.' + os.path.basename(input_file) + '/' + os.path.basename(input_file)
    return output_file_name

def return_primary_command(args, output_file_name, input_file, input_json):
    primary_command = 'java -Dconfig.file=' + re.sub('input.json', 'overrides.conf', input_json) + ' -jar ' + args.cromwell_path + ' run ' + re.sub('input.json', 'haplotypecaller-gvcf-gatk4.wdl', input_json) + ' -i ' + input_json
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
    args = parse_args()
    clean_arg_paths(args)
    
    os.system('module load gcc/6.2.0 python/3.6.0 java')  

    input_files = return_input_files(args, 'bam') if args.input_directory is not None else [args.input_file_path]
    input_files = sort_by_size(input_files)
    
    for input_file in input_files:

        input_json = generate_input_json(args, input_file)
        slurm_command = return_slurm_command(args)
        output_file_name = gen_output_file_name(args, input_file)
        primary_command = return_primary_command(args, output_file_name, input_file, input_json)

        sh_file_name = gen_sh_file_name(args, output_file_name)
        write_out(args, slurm_command, primary_command, sh_file_name)

        sample_name = ntpath.basename(sh_file_name).split('.bam')[0] + '.g.vcf.gz'  
        path_to_vcf = os.path.join(ntpath.dirname(ntpath.dirname(ntpath.dirname(sh_file_name))), sample_name)

        if not os.path.isfile(path_to_vcf):
            #print(path_to_vcf)
            submit_job(sh_file_name)

        
if __name__ == "__main__":
    main()
