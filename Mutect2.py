

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
from time import sleep

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
    parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf', help='path to dbsnp file')
    # parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_hg38_common_all_20160601.vcf', help='path to dbsnp file')
    parser.add_argument('-gnomad', '--gnomad_path', default='/n/data1/hms/dbmi/park/victor/software/GATK_bundle/af-only-gnomad.hg19.vcf', help='path to cosmic file' )
    parser.add_argument('-scatter', '--scatter_size', default='50')
    return parser.parse_args()

def clean_arg_paths(args):
    """Modifies all user-inputted directory paths such that they end with a '/'"""
    d = vars(args)
    for arg in d.keys():
        if 'output_directory' in arg and d[arg]=='./': d[arg] = os.getcwd() + '/.' + ntpath.basename(args.input_tumor_path) + '_v_' + ntpath.basename(args.input_normal_path)
        if 'directory' in arg and d[arg] is not None and d[arg][-1] != '/': d[arg] += '/'

def generate_regions_files(args):
    os.makedirs(os.path.dirname(args.output_directory + '.Mutect2/.regions/'), exist_ok=True)
    os.system(args.gatk4_path + ' SplitIntervals' + '\\' + '\n' + \
     '\t' + '-R ' + args.reference_path + ' \\' + '\n' + \
     '\t' + '-scatter ' + args.scatter_size + ' \\' + '\n' + \
     '\t' + '-O ' + args.output_directory + '.Mutect2/.regions/' + ' \\')
    sleep(15)

def return_region_files(args):
    region_files = [args.output_directory + '.Mutect2/.regions/' + file for file in os.listdir(args.output_directory + '.Mutect2/.regions/')]
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

    generate_regions_files(args)
    region_files = return_region_files(args)

    for region_file in region_files:
        slurm_command = return_slurm_command(args)
        output_file_name = gen_output_file_name(args, region_file)

        if args.panel == 'nopath':
            primary_command = return_new_primary_command(args, output_file_name, region_file)
        else: 
            primary_command = return_pon_command(args, output_file_name, region_file)

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
            submit_job(sh_file_name)
        #submit_job(sh_file_name)
"""
        sample_name = ntpath.basename(sh_file_name).split('.bam')[0]
        path_to_vcf = os.path.join('/n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK/tumor_bams/.Mutect2', sample_name)
        vcf = ntpath.basename(sh_file_name).split('.sh')[0]
        path_to_vcf = os.path.join(path_to_vcf, vcf)
        if not os.path.isfile(path_to_vcf):
            #print(sh_file_name)
            submit_job(sh_file_name)
"""
if __name__ == "__main__":
    main()
