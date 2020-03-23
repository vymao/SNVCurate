"""This script takes in text files of lists of BAM files and executes the GATK Haplotype and Joint Mutation Callers to produce the desired Output."""

"""
python3 /n/data1/hms/dbmi/park/victor/scripts/other/MuSE_run.py -tumor /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned -normal /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/MuSE_v2 -csv /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned/PancSeq_WES_cohort.csv -mode call -r2 20 -data_type WES

python3 /n/data1/hms/dbmi/park/victor/scripts/other/MuSE_run.py -tumor /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/MuSE_v2 -csv /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/MuSE_v2/PancSeq_WES_cohort/.MuSE/called_list.csv -mode sump -data_type WES -t 0-01:00:00 --mem_per_cpu 5G -r1 0 -r2 1


python3 /n/data1/hms/dbmi/park/victor/scripts/other/MuSE_run.py -tumor /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel -normal /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel/MuSE -csv /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel/Panc_new_matched_samples_WES.csv -mode call -r1 2 
python3 /n/data1/hms/dbmi/park/victor/scripts/other/MuSE_run.py -tumor /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel/MuSE/Panc_new_matched_samples_WES -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel/MuSE -csv /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel/MuSE/Panc_new_matched_samples_WES/tumor_list.csv -mode sump -t 0-01:00:00 --mem_per_cpu 5G -r1 0 -r2 1



"""
import argparse
import os
import re
import ntpath
import sys

def parse_args():
    """Uses argparse to enable user to customize script functionality""" 
    #BE SURE TO ADD MULTIPLE TEXT FILE ARGUMENTS Later
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-tumor', '--input_tumor_path', help='path to input tumor file')
    parser.add_argument('-normal', '--input_normal_path', help='path to normal file')
    parser.add_argument('-csv', help='csv containing matched tumor/normal pairs')
    parser.add_argument('-out', '--output_directory', default='./', help='directory to which the output directory "/.Mutect2/" will be written to')
    parser.add_argument('-n', '--num_cores', default='1', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='2-00:00:00', help='slurm job submission option')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='15G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', default='victor_mao@hms.harvard.edu', help='slurm job submission option')
    parser.add_argument('-reference', '--reference_path', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta', help='path to reference_path file')
    # parser.add_argument('-reference', '--reference_path', default='/n/dlsata1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38/Homo_sapiens_assembly38.fasta', help='path to reference_path 
    parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf.gz', help='path to dbsnp file')
    # parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_hg38_common_all_20160601.vcf', help='path to dbsnp file')
    parser.add_argument('-mode', default='call')
    parser.add_argument('-data_type', default='WGS')
    parser.add_argument('-r1', default=1, help='Lower range bound of indices of BAMs to run')
    parser.add_argument('-r2', default=100000, help='Upper range bound of indices of BAMs to run')

    return parser.parse_args()

def arg_clean(args):
    """Cleaning the parsed arguments. Specifically, this function does the following: 
    1. Finds the host filename
    2. Finds the desired output directory
    3. Creates a new output directory within the specified directory
    4. Specifies the appropriate script, depending on the desired pipeline
    """
    if '.csv' in args['csv']:
        filename = re.findall('/[A-Za-z0-9_]*\.', args['csv'])[0][1:-1]
    else:
        args['csv'] = args['csv'] + '.csv'
        filename = re.findall('/[A-Za-z0-9_]*\.', args['csv'])[0][1:-1]

    #if args['out'] == './': output_dir = args['out']
    if args['output_directory'] == '.' or args['output_directory'] == './': output_dir = os.getcwd() 
    if args['mode'] == 'call':       
        output_dir = os.path.join(args['output_directory'], filename)
        os.makedirs(output_dir, exist_ok = True)
    else: 
        output_dir = args['output_directory']
    print(output_dir)
    return output_dir

#This needs further editing to accomodate more relevant scripts
"""
    if args['pipeline'].lower() == 'somatic':
        pipeline = 'GATK_Somatic_SNPs_Indels/Mutect2.py'
    else:
        if args['queue'] == 'priopark':
            pipeline = 'GATK_Germline_SNPs_Indels/HaplotypeCaller_original.py'
        else:
            pipeline = 'GATK_Germline_SNPs_Indels/HaplotypeCaller.py'

    return output_dir, pipeline
"""
def get_column(csv, sample):
    with open(csv, 'r') as file:
        for index, line in enumerate(file):
            if index == 0:
                result = [x.strip().lower() for x in line.split(',')]
                if sample == "T":
                    for header in result: 
                        if "_t" in header or "tum" in header:
                            return result.index(header)
                else:
                    for header in result: 
                        if "_n" in header or "norm" in header:
                            return result.index(header)

def get_bam(csv, row, column):
    with open(csv, 'r') as file:
        for idx, line in enumerate(file):
            if idx == row: 
                result = [x.strip() for x in line.split(',')]
                if result[column] != 'N/A':
                    return result[column]

def main(): 
    tools_dir = '/n/data1/hms/dbmi/park/victor/scripts/other/MuSE.py'
    args = vars(parse_args())
    output_dir = arg_clean(args)
    

    if int(args['r1']) == 0 and args['mode'].lower() == 'call':
        print("Invalid r1 parameter. r1 must be greater than 0 to account for headers")
        sys.exit()

    os.system('module load gcc/6.2.0 python/3.6.0 java perl')

    if args['mode'].lower() == 'call':
        tumor_index = get_column(args['csv'], "T")
        normal_index = get_column(args['csv'], "N")
        with open(args['csv'], 'r') as f:
            for index, line in enumerate(f):
                if not line.isspace() and index in range(int(args['r1']), int(args['r2'])):
                    current_tumor_sample = get_bam(args['csv'], index, tumor_index)
                    current_normal_sample = get_bam(args['csv'], index, normal_index)
                    #print(current_normal_sample)

                    if current_tumor_sample is None or current_normal_sample is None:
                        continue
                    else: 
                        tumor_sample = os.path.join(args['input_tumor_path'], get_bam(args['csv'], index, tumor_index))
                        normal_sample = os.path.join(args['input_normal_path'], get_bam(args['csv'], index, normal_index))

                        os.system('python3 ' + tools_dir + ' -tumor ' + tumor_sample + ' -normal ' + normal_sample + ' -out ' + output_dir + ' -t ' + args['runtime'] + ' -n ' + args['num_cores'] + 
                            ' -p ' + args['queue'] + ' --mail_user ' + args['mail_user'] + ' --mem_per_cpu ' + args['mem_per_cpu'] + ' --mail_type ' + args['mail_type'] + ' -reference ' + args['reference_path'] + ' -mode ' + args['mode'] 
                            + ' -data_type ' + args['data_type'])

                        #print('python3 ' + tools_dir + ' -tumor ' + tumor_sample + ' -normal ' + normal_sample + ' -out ' + output_dir + ' -t ' + args['runtime'] + 
                        #    ' -p ' + args['queue'] + ' --mail_user ' + args['mem_per_cpu'] + ' -reference ' + args['reference_path'] + ' -mode ' + args['mode'])

    if args['mode'].lower() == 'sump':
        with open(args['csv'], 'r') as f:
            for index, line in enumerate(f):
                if not line.isspace() and index in range(int(args['r1']), int(args['r2'])):
                    tumor = line.strip();
                    os.system('python3 ' + tools_dir + ' -tumor ' + tumor + ' -out ' + output_dir + ' -t ' + args['runtime'] + ' -n ' + args['num_cores'] + 
                            ' -p ' + args['queue'] + ' --mail_user ' + args['mail_user'] + ' --mem_per_cpu ' + args['mem_per_cpu'] + ' --mail_type ' + args['mail_type'] + ' -reference ' + args['reference_path'] + ' -mode ' + args['mode'] 
                            + ' -data_type ' + args['data_type'])
if __name__ == "__main__":
    main()

