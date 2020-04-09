

"""
python3 /n/data1/hms/dbmi/park/victor/scripts/other/MuSE.py -tumor -normal -out -mode call 

python3 /n/data1/hms/dbmi/park/victor/scripts/other/MuSE.py -tumor /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/edited_bams/QRF116304_reordered_contigs.bam -normal /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/edited_bams/normals/wgs_normal.bam -out /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/MuSE -reference /home/mk446/BiO/Install/GATK-bundle/2.8/hg19/ucsc.hg19.fasta -mode call 


python3 /n/data1/hms/dbmi/park/victor/scripts/other/MuSE.py -tumor /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/MuSE/.MuSE/QRF116304_reordered_contigs.bam_v_wgs_normal.bam.MuSE.txt -out /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/MuSE -reference /home/mk446/BiO/Install/GATK-bundle/2.8/hg19/ucsc.hg19.fasta -mode sump 

python3 /n/data1/hms/dbmi/park/victor/scripts/other/MuSE.py -tumor /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/edited_bams/QRF116304_reordered_contigs.bam -normal /n/data1/hms/dbmi/park/vinay/analysis/externalDatasets/SRP145073_lee_etal_nature_2018_glioblastoma_svz_spatial/alignment/20191121_readgroup/SRR7138430.sorted.readgroup.bam -out /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/MuSE -reference /home/mk446/BiO/Install/GATK-bundle/2.8/hg19/ucsc.hg19.fasta -mode call -data_type WES

python3 /n/data1/hms/dbmi/park/victor/scripts/other/MuSE.py -tumor /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/MuSE/.MuSE/QRF116304_reordered_contigs.bam_v_SRR7138430.sorted.readgroup.bam.MuSE.txt -out /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/MuSE -reference /home/mk446/BiO/Install/GATK-bundle/2.8/hg19/ucsc.hg19.fasta -mode sump -t 0-02:00:00 --mem_per_cpu 1G

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
    parser.add_argument('-tumor', '--input_tumor_path', help='path to input tumor file')
    parser.add_argument('-normal', '--input_normal_path', help='path to normal file')
    parser.add_argument('-out', '--output_directory', default='./', help='directory to which the output directory "/.MuSE/" will be written to')
    parser.add_argument('-n', '--num_cores', default='1', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='1-00:00:00', help='slurm job submission option')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='5G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', default='victor_mao@hms.harvard.edu', help='slurm job submission option')
    parser.add_argument('-reference', '--reference_path', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta', help='path to reference_path file')
    # parser.add_argument('-reference', '--reference_path', default='/n/dlsata1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38/Homo_sapiens_assembly38.fasta', help='path to reference_path 
    parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf.gz', help='path to dbsnp file')
    #parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_hg19_common_all_20160601.vcf.gz', help='path to dbsnp file')
    parser.add_argument('-mode')
    parser.add_argument('-data_type', default='WGS')

    parser.add_argument('-cn', default="1", help='number of cores for Cromwell jobs')
    parser.add_argument('-ct', default="1000", help='cromwell run time; please specify as number of minutes')
    parser.add_argument('-cm', default='7000', help='cromwell cpu memory per core')
    parser.add_argument('-cromwell', '--cromwell_path', default='/n/shared_db/singularity/hmsrc-gatk/cromwell-43.jar', help='path to cromwell.jar file')
    return parser.parse_args()

def main():
    args = parse_args()
    clean_arg_paths(args)

    dirname = os.path.dirname(os.path.abspath(__file__))
    overrides = os.path.join(dirname, 'Overrides.config')
    wdl = os.path.join(dirname, 'MuSE.wdl')
    json = os.path.join(dirname, 'MuSE_Inputs.json') 

    input_json, input_config, input_wdl = generate_cromwell_inputs(args, json, wdl, overrides)

    slurm_command = return_slurm_command(args)
    output_file_name = gen_output_file_name(args)

    primary_command = return_primary_command(args, output_file_name, input_json, input_config, input_wdl)

    sh_file_name = gen_sh_file_name(args, output_file_name)
    write_out(args, slurm_command, primary_command, sh_file_name)
    sample_name = ntpath.basename(sh_file_name).split('.bam')[0]
    muse_path = ntpath.dirname(ntpath.dirname(sh_file_name))
    #path_to_vcf = os.path.join(mutect_path, sample_name)
    path_to_vcf = muse_path
    vcf = ntpath.basename(sh_file_name).split('.sh')[0]
    dirname = vcf.split('.')[0]
    path_to_vcf = os.path.join(path_to_vcf, os.path.join(dirname, vcf)) 
    #path_to_vcf = os.path.join(path_to_vcf, vcf)
    #print(path_to_vcf)
    if not os.path.isfile(path_to_vcf):
        print(path_to_vcf)
        #print(primary_command)
        #submit_job(sh_file_name)

    """
    if args.mode == 'call':
        primary_command = return_primary_call_command(args, output_file_name)
    else:
        primary_command = return_primary_sump_command(args, output_file_name)

    sh_file_name = gen_sh_file_name(args, output_file_name)
    write_out(args, slurm_command, primary_command, sh_file_name)

    sample_name = ntpath.basename(sh_file_name).split('.bam')[0]
    mutect_path = ntpath.dirname(ntpath.dirname(sh_file_name))
    #path_to_vcf = os.path.join(mutect_path, sample_name)
    path_to_vcf = mutect_path
    vcf = ntpath.basename(sh_file_name).split('.sh')[0]
    dirname = vcf.split('.')[0]
    #path_to_vcf = os.path.join(path_to_vcf, os.path.join(dirname, vcf)) 
    path_to_vcf = os.path.join(path_to_vcf, vcf)
    #print(path_to_vcf)
    if not os.path.isfile(path_to_vcf):
        #print(path_to_vcf)
        #print(primary_command)
        submit_job(sh_file_name)
    """

def clean_arg_paths(args):
    """Modifies all user-inputted directory paths such that they end with a '/'"""
    d = vars(args)
    for arg in d.keys():
        if 'output_directory' in arg and d[arg]=='./': d[arg] = os.getcwd() + '/.' + ntpath.basename(args.input_tumor_path) + '_v_' + ntpath.basename(args.input_normal_path)
        if 'directory' in arg and d[arg] is not None and d[arg][-1] != '/': d[arg] += '/'


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

def gen_output_file_name(args):
    sample = ntpath.basename(args.input_tumor_path).split('.')[0]
    output_file_name = args.output_directory + '.MuSE/' + sample + '.vcf'
    """
    if args.mode == 'call':
        output_file_name = args.output_directory + '.MuSE/' + ntpath.basename(args.input_tumor_path) + '_v_' + ntpath.basename(args.input_normal_path)
    else: 
        sample = ntpath.basename(args.input_tumor_path).split('.')[0]
        output_file_name = args.output_directory + '.MuSE/' + sample + '.vcf'
    """
    return output_file_name

def generate_cromwell_inputs(args, json_file, wdl, overrides):
    input_file = args.input_tumor_path
    dir = args.output_directory + '.MuSE/' + '.' + os.path.basename(input_file) + '/'
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
        d["MuSE.input_bam"] = input_file
        d["MuSE.input_bam_index"] = path
        d["MuSE.output_directory"] = os.path.join(args.output_directory,'.MuSE/' + bam_sample.split('.')[0] + '/')
        d["MuSE.ref_dict"] = os.path.join(dict_path, ref + '.dict')
        d["MuSE.ref_fasta"] = args.r
        d["MuSE.ref_fasta_index"] = args.r + '.fai'
        d["MuSE.gatk_path"] = args.gatk
        d["MuSE.dbSNP"] = args.dbsnp_path
        d["MuSE.data_type"] = args.data_type

   
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

def return_primary_command(args, output_file_name, input_json, input_config, input_wdl):
    primary_command = 'java -Dconfig.file=' + input_config + ' -jar ' + args.cromwell_path + ' run ' + input_wdl + ' -i ' + input_json
    return primary_command

"""
def return_primary_call_command(args, output_file_name):
    primary_command = 'MuSE call' + \
    ' -f ' + args.reference_path + \
    ' -O ' + output_file_name + \
    ' ' + args.input_tumor_path + \
    ' ' + args.input_normal_path
    return primary_command

def return_primary_sump_command(args, output_file_name):
    primary_command = 'MuSE sump' + \
    ' -I ' + args.input_tumor_path + \
    ' -O ' + output_file_name + \
    ' -D ' + args.dbsnp_path
    if args.data_type == 'WGS': 
        primary_command += ' -G'
    else:
        primary_command += ' -E'
    return primary_command
"""
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
