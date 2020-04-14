

"""
python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2.py -tumor  /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK/s_DS_bkm_129_T1_bc0028_Proj_5065_E_L000_mrg_cl_aln_srt_MD_IR_BR.bam -normal /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK/s_DS_bkm_129_N_bc0072_Proj_5065_E_L000_mrg_cl_aln_srt_MD_IR_BR.bam -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/Mutect_recalibrated 

python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2.py -tumor  /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK/s_DS_bkm_129_T2_bc0029_Proj_5065_E_L000_mrg_cl_aln_srt_MD_IR_BR.bam -normal /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK/s_DS_bkm_129_N_bc0072_Proj_5065_E_L000_mrg_cl_aln_srt_MD_IR_BR.bam -out /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK/tumor_bams_cp


python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2.py -tumor /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned/14408_CCPM_0400245_14408_CCPM_0400245_T1_SM-DFY6G_RECAL.bam -normal /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned/14408_CCPM_0400245_14408_CCPM_0400245_BLOOD_SM-DRXSH_RECAL.bam -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/reference_test/PancSeq_WES_cohort -t 2-0:00:00 -p medium -reference /n/data1/hms/dbmi/park/victor/references/human_g1k_v37.fasta


python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2.py -tumor /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned/400069_0400069_T1_SM-AJDYG_RECAL.bam -normal /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned/400069_0400069_N_SM-AJDZA_RECAL.bam -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/reference_test/PancSeq_WES_cohort -t 2-0:00:00 -p medium -reference /n/data1/hms/dbmi/park/victor/references/human_g1k_v37.fasta

python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2.py -tumor /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/edited_bams/QRF116304.reordered_contigs.bam -normal /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned/400075_400075_N_2_SM-AX6C9_RECAL.bam -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/reference_test/PancSeq_WES_cohort -t 2-0:00:00 -p medium -reference /n/data1/hms/dbmi/park/victor/references/human_g1k_v37.fasta


"""


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
    parser = argparse.ArgumentParser()
    parser.add_argument('-tumor', '--input_tumor_path', help='path to input tumor file')
    parser.add_argument('-normal', '--input_normal_path', help='path to normal file')
    parser.add_argument('-pon', '--panel', default='nopath', help='path to panel of norms')
    parser.add_argument('-out', '--output_directory', default='./', help='directory to which the output directory "/.Mutect2/" will be written to')
    parser.add_argument('-n', '--num_cores', default='2', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='3-00:00:00', help='slurm job submission option')
    parser.add_argument('-p', '--queue', default='medium', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='20G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', default='victor_mao@hms.harvard.edu', help='slurm job submission option')
    parser.add_argument('-gatk', '--gatk_path', default='/home/mk446/BiO/Install/GATK3.5_160425_g7a7b7cd/GenomeAnalysisTK.jar', help='path to software')
    parser.add_argument('-gatk_new', '--gatk_path_new', default='/home/mk446/BiO/Install/GATK4.1.2.0/gatk', help='path to software')
    parser.add_argument('-gatk4', '--gatk4_path', default='/n/data1/hms/dbmi/park/alon/software/gatk/gatk-4.0.3.0/gatk', help='path to software')
    parser.add_argument('-reference', '--reference_path', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta', help='path to reference_path file')
    # parser.add_argument('-reference', '--reference_path', default='/n/dlsata1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38/Homo_sapiens_assembly38.fasta', help='path to reference_path 
    parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/af-only-gnomad.raw.sites.b37.vcf', help='path to dbsnp file')
    #parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_hg38_common_all_20160601.vcf', help='path to dbsnp file')
    parser.add_argument('-gnomad', '--gnomad_path', default='/n/data1/hms/dbmi/park/victor/software/GATK_bundle/af-only-gnomad.raw.sites.b37.vcf.gz', help='path to cosmic file' )
    #parser.add_argument('-scatter', '--scatter_size', default='50')
    parser.add_argument('-interval_list', default='/n/data1/hms/dbmi/park/victor/software/MuTecT2_b37_scattered_intervals.txt')

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
    wdl = os.path.join(dirname, 'MuTecT.wdl')
    json = os.path.join(dirname, 'MuTecT_Inputs.json') 

    sample_name = ntpath.basename(args.input_tumor_path).split('.bam')[0]
    path_to_vcf = os.path.join(os.path.join(args.output_directory, os.path.join(".MuTecT2", sample_name)), sample_name + '.vcf') 
    vcf_dir = os.path.join(args.output_directory, os.path.join(".MuTecT2", sample_name))
    os.makedirs(vcf_dir, exist_ok=True) 

    slurm_command = return_slurm_command(args)
    output_file_name = gen_output_file_name(args)

    """
    if args.panel == 'nopath':
        primary_command = return_new_primary_command(args, output_file_name, region_file)
    else: 
        primary_command = return_pon_command(args, output_file_name, region_file)
    """
    input_json, input_config, input_wdl = generate_cromwell_inputs(args, json, wdl, overrides)
    primary_command = return_primary_command(args, output_file_name, input_json, input_config, input_wdl)
    
    sh_file_name = gen_sh_file_name(args, output_file_name)
    write_out(args, slurm_command, primary_command, sh_file_name)

    if not os.path.isfile(path_to_vcf):
        #print(path_to_vcf)
        submit_job(sh_file_name)
    #submit_job(sh_file_name)

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
    slurm_command += 'module load gcc/6.2.0 java/jdk-1.8u112 bcftools/1.9' + '\n'
    return slurm_command

def generate_cromwell_inputs(args, json_file, wdl, overrides):
    input_file = args.input_tumor_path
    normal_file = args.input_normal_path
    dir = args.output_directory + '.MuTecT2/' + '.' + os.path.basename(input_file).split('.')[0] + '/'
    os.makedirs(dir, exist_ok=True)

    bam_dir = os.path.dirname(input_file)
    bam_sample = os.path.basename(input_file)
    normal_dir = os.path.dirname(normal_file)
    normal_sample = os.path.basename(normal_file)

    bai_suffix = '.bai'
    bam_path = os.path.join(bam_dir, re.sub('.bam', '.bam.bai', bam_sample))
    normal_path = os.path.join(normal_dir, re.sub('.bam', '.bam.bai', normal_sample))

    if not (os.path.isfile(bam_path) and os.access(bam_path, os.R_OK)):
        bam_path = os.path.join(bam_dir, re.sub('.bam', '.bai', bam_sample))
    if not (os.path.isfile(normal_path) and os.access(normal_path, os.R_OK)):
        normal_path = os.path.join(normal_dir, re.sub('.bam', '.bai', normal_sample))
    
    copyfile(json_file, dir + 'Input.json')

    dict_path = os.path.dirname(args.reference_path)
    ref = os.path.basename(args.reference_path).split('.fa')[0]
    
    with open(dir + 'Input.json') as f:
        data = f.read()
        d = json.loads(data)
        d["MuTecT.input_bam"] = input_file
        d["MuTecT.input_bam_index"] = bam_path
        d["MuTecT.output_directory"] = os.path.join(args.output_directory,'.MuTecT2/' + bam_sample.split('.')[0] + '/')
        d["MuTecT.ref_dict"] = os.path.join(dict_path, ref + '.dict')
        d["MuTecT.ref_fasta"] = args.reference_path
        d["MuTecT.ref_fasta_index"] = args.reference_path + '.fai'
        d["MuTecT.gnomad"] = args.gnomad_path
        d["MuTecT.gatk_path"] = args.gatk_path_new
        d["MuTecT.interval_list"] = args.interval_list
        d["MuTecT.gnomad_index"] = args.gnomad_path + '.tbi'

        if args.panel == "nopath":
            d["MuTecT.normal_bam"] = normal_file
            d["MuTecT.normal_bam_index"] = normal_path
            d["MuTecT.mode"] = "normal"
        else: 
            d["MuTecT.panel"] = args.panel
            d["MuTecT.panel_bam_index"] = args.panel + '.idx'
            d["MuTecT.mode"] = "panel"

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

def gen_output_file_name(args):
    sample = ntpath.basename(args.input_tumor_path).split('.')[0]
    output_file_name = args.output_directory + '.MuTecT2/.' + os.path.basename(args.input_tumor_path).split('.')[0] + '/' + sample + '.vcf'

    return output_file_name

"""
def gen_output_file_name(args, region_file):
    if args.panel == 'nopath':
        output_file_name = args.output_directory + '.Mutect2/' + ntpath.basename(args.input_tumor_path) + '_v_' + ntpath.basename(args.input_normal_path) + '_' + ntpath.basename(region_file[1:]) + '.vcf'
    else: 
        output_file_name = args.output_directory + '.Mutect2/' + ntpath.basename(args.input_tumor_path) + '_v_' + ntpath.basename(args.panel) + '_' + ntpath.basename(region_file[1:]) + '.vcf'

    return output_file_name

def return_primary_command(args, output_file_name, region_file):
    primary_command = 'java -jar ' + args.gatk_path + ' \\' + '\n' + \
     '\t' + '-T MuTect2' + ' \\' + '\n' + \
     '\t' + '-R ' + args.reference_path + ' \\' + '\n' + \
     '\t' + '-I:tumor ' + args.input_tumor_path  + ' \\' + '\n' + \
     '\t' + '-I:normal ' + args.input_normal_path + ' \\' + '\n' + \
     '\t' + '--dbsnp ' + args.dbsnp_path + ' \\' + '\n' + \
     '\t' + '-L ' + args.output_directory + '.Mutect2/.regions/' + ntpath.basename(region_file) + ' \\' + '\n' + \
     '\t' + '-dt NONE' + ' \\' + '\n' + \
     '\t' + '-ploidy 2' + ' \\' + '\n' + \
     '\t' + '-o ' + output_file_name
    return primary_command

def return_new_primary_command(args, output_file_name, region_file):
    primary_command = args.gatk_path_new + ' Mutect2 ' + '\\' + '\n' + \
     '\t' + '-R ' + args.reference_path + ' \\' + '\n' + \
     '\t' + '-I ' + args.input_tumor_path  + ' \\' + '\n' + \
     '\t' + '-normal ' + args.input_normal_path + ' \\' + '\n' + \
     '\t' + '--germline-resource ' + args.gnomad_path + ' \\' + '\n' + \
     '\t' + '-L ' + args.output_directory + '.Mutect2/.regions/' + ntpath.basename(region_file) + ' \\' + '\n' + \
     '\t' + '-O ' + output_file_name
    return primary_command

def return_pon_command(args, output_file_name, region_file):
    primary_command = 'java -jar ' + args.gatk_path + ' \\' + '\n' + \
     '\t' + '-T MuTect2' + ' \\' + '\n' + \
     '\t' + '-R ' + args.reference_path + ' \\' + '\n' + \
     '\t' + '-I:tumor ' + args.input_tumor_path  + ' \\' + '\n' + \
     '\t' + '-PON ' + args.panel + ' \\' + '\n' + \
     '\t' + '--dbsnp ' + args.dbsnp_path + ' \\' + '\n' + \
     '\t' + '-L ' + args.output_directory + '.Mutect2/.regions/' + ntpath.basename(region_file) + ' \\' + '\n' + \
     '\t' + '-dt NONE' + ' \\' + '\n' + \
     '\t' + '-ploidy 2' + ' \\' + '\n' + \
     '\t' + '-o ' + output_file_name
    return primary_command

def return_new_pon_command(args, output_file_name, region_file):
    primary_command = args.gatk_path_new + ' Mutect2 ' + '\\' + '\n' + \
     '\t' + '-R ' + args.reference_path + ' \\' + '\n' + \
     '\t' + '-I ' + args.input_tumor_path  + ' \\' + '\n' + \
     '\t' + '--panel-of-normals ' + args.panel + ' \\' + '\n' + \
     '\t' + '--germline-resource ' + args.gnomad_path + ' \\' + '\n' + \
     '\t' + '-L ' + args.output_directory + '.Mutect2/.regions/' + ntpath.basename(region_file) + ' \\' + '\n' + \
     '\t' + '-O ' + output_file_name
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
