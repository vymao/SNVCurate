from time import sleep
import argparse
import os
import re
import ntpath
import glob
import sys
import fileinput
from shutil import copyfile
import json


def parse_args():
    """Uses argparse to enable user to customize script functionality"""
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in_dir', '--input_directory', help='path to directory containing input files')
    parser.add_argument('-in_file', '--input_file_path', help='path to input file')
    parser.add_argument('-out', '--output_directory', default='./', help='directory to which the output directory "/.HaplotypeCaller/" will be written to')
    parser.add_argument('-n', '--num_cores', default='1', help='slurm job submission option')
    parser.add_argument('-cn', default="1", help='number of cores for Cromwell jobs')
    parser.add_argument('-t', '--runtime', default='2-00:00:00', help='slurm job submission option')
    parser.add_argument('-ct', default="3000", help='cromwell run time; please specify as number of minutes')
    parser.add_argument('-p', '--queue', default='long', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='15G', help='slurm job submission option')
    parser.add_argument('-cm', default='5000', help='cromwell cpu memory per core')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', default='victor_mao@hms.harvard.edu', help='slurm job submission option')
    parser.add_argument('-cromwell', '--cromwell_path', default='/n/shared_db/singularity/hmsrc-gatk/cromwell-43.jar', help='path to cromwell.jar file')
    parser.add_argument('-r', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta')
    parser.add_argument('-gatk', default='/home/mk446/BiO/Install/GATK4.1.2.0//gatk', help='path to software')
    parser.add_argument('-reference_name', default='b37', help='hg19, b37, etc.')
    parser.add_argument('-picard', default='/home/mk446/BiO/Install/picard-tools-2.5.0/picard.jar')
    parser.add_argument('-read_groups', default='single')
    return parser.parse_args()

def main():
    args = parse_args()
    clean_arg_paths(args)

    dirname = os.path.dirname(os.path.abspath(__file__))
    overrides = os.path.join(dirname, 'Overrides.config')
    wdl = os.path.join(dirname, 'Germline_workflow.wdl')
    json = os.path.join(dirname, 'Germline_workflow.json') 

    input_files = return_input_files(args, 'bam') if args.input_directory is not None else [args.input_file_path]
    input_files = sort_by_size(input_files)
    
    for input_file in input_files:
        sample = os.path.basename(input_file).split('.')[0]
        sample_dir = os.path.join(args.output_directory,'.HaplotypeCaller/' + sample + '/')
        os.makedirs(sample_dir, exist_ok=True)

        input_json, input_config, input_wdl = generate_cromwell_inputs(args, input_file, json, wdl, overrides)
        slurm_command = return_slurm_command(args)
        output_file_name = gen_output_file_name(args, input_file)
        primary_command = return_primary_command(args, output_file_name, input_file, input_json, input_config, input_wdl)

        sh_file_name = gen_sh_file_name(args, output_file_name)
        write_out(args, slurm_command, primary_command, sh_file_name)

        sample_name = ntpath.basename(sh_file_name).split('.bam')[0] + '.vcf'  
        path_to_vcf = os.path.join(ntpath.dirname(ntpath.dirname(ntpath.dirname(sh_file_name))), sample_name)

        if not os.path.isfile(path_to_vcf):
            submit_job(sh_file_name)


def clean_arg_paths(args):
    """Modifies all user-inputted directory paths such that they end with a '/'"""
    d = vars(args)
    for arg in d.keys():
        if 'input_directory' in arg and d[arg]=='./': d[arg] = os.getcwd()
        if 'output_directory' in arg and d[arg]=='./': d[arg] = os.getcwd()    
        if 'directory' in arg and d[arg] is not None and d[arg][-1] != '/': d[arg] += '/'
    output_dir = re.sub(" ", "", d["output_directory"])
    if output_dir[len(output_dir) - 1] is not "/":
        d["output_directory"] = output_dir + "/"

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

def generate_cromwell_inputs(args, input_file, json_file, wdl, overrides):
    dir = args.output_directory + '.HaplotypeCaller/' + '.' + os.path.basename(input_file).split('.')[0] + '/'
    os.makedirs(dir, exist_ok=True)

    bam_dir = os.path.dirname(input_file)
    bam_sample = os.path.basename(input_file)

    bai_suffix = '.bai'
    path = os.path.join(bam_dir, re.sub('.bam', '.bam.bai', bam_sample))
    if not (os.path.isfile(path) and os.access(path, os.R_OK)):
        path = os.path.join(bam_dir, re.sub('.bam', '.bai', bam_sample))
    
    copyfile(json_file, dir + 'Input.json')

    dict_path = os.path.dirname(args.r)
    ref = os.path.basename(args.r).split('.fa')[0]
    
    with open(dir + 'Input.json') as f:
        data = f.read()
        d = json.loads(data)
        d["HaplotypeCallerGvcf_GATK4.input_bam"] = input_file
        d["HaplotypeCallerGvcf_GATK4.input_bam_index"] = path
        d["HaplotypeCallerGvcf_GATK4.output_directory"] = os.path.join(args.output_directory,'.HaplotypeCaller/' + bam_sample.split('.')[0] + '/')
        d["HaplotypeCallerGvcf_GATK4.ref_dict"] = os.path.join(dict_path, ref + '.dict')
        d["HaplotypeCallerGvcf_GATK4.ref_fasta"] = args.r
        d["HaplotypeCallerGvcf_GATK4.ref_fasta_index"] = args.r + '.fai'
        d["HaplotypeCallerGvcf_GATK4.gatk_path"] = args.gatk
        d["HaplotypeCallerGvcf_GATK4.picard_path"] = args.picard
        d["HaplotypeCallerGvcf_GATK4.bam_directory"] = os.path.dirname(input_file)
        
        if args.read_groups == "single":
            d["HaplotypeCallerGvcf_GATK4.read_groups"] = "single"
        elif args.read_groups == "multiple":
            d["HaplotypeCallerGvcf_GATK4.read_groups"] = "multiple"

        if args.reference_name.lower() is not "b37":
            d["HaplotypeCallerGvcf_GATK4.reference"] = args.reference_name

   
    with open(dir + 'Input.json', 'w') as f:
        f.write(json.dumps(d))

    copyfile(overrides, dir + 'Overrides.config')

    with fileinput.FileInput(dir + 'Overrides.config', inplace=True) as file:
        for line in file:
            print(line.replace(
            "medium", args.queue).replace("!@#$", args.ct).replace("%^&*", args.cm).replace("Int cpus = 1", "Int cpus = " + args.cn), end='')
        

    if args.queue == 'park' or args.queue == 'priopark': 
        f = open(dir + 'Overrides.config', "r")
        contents = f.readlines()
        f.close()

        contents.insert(95, "            --account=${account_name} \\\n")

        f = open(dir + 'Overrides.config', "w")
        contents = "".join(contents)
        f.write(contents)
        f.close()

    copyfile(wdl, os.path.join(dir, 'workflow.wdl'))

    return os.path.join(dir,'Input.json'), os.path.join(dir, 'Overrides.config'), os.path.join(dir, 'workflow.wdl')

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
    slurm_command += 'module load gcc/6.2.0 java/jdk-1.8u112 bcftools/1.9' + '\n'
    return slurm_command

def gen_output_file_name(args, input_file):
    sample = os.path.basename(input_file).split('.')[0]
    output_file_name = args.output_directory + '.HaplotypeCaller/' + '.' + sample + '/' + sample
    return output_file_name

def return_primary_command(args, output_file_name, input_file, input_json, input_config, input_wdl):
    primary_command = 'java -Dconfig.file=' + input_config + ' -jar ' + args.cromwell_path + ' run ' + input_wdl + ' -i ' + input_json
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

        
if __name__ == "__main__":
    main()
