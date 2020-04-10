"""This script takes in text files of lists of BAM files and executes the GATK Haplotype and Joint Mutation Callers to produce the desired Output."""

"""
Testing queries:
python3 /n/data1/hms/dbmi/park/victor/scripts/BAMs_txt.py -input_path /n/data1/hms/dbmi/park/victor/other/random/Testfile.txt -output_path /n/data1/hms/dbmi/park/victor/other/tests/ -pipeline Germline
python3 /n/data1/hms/dbmi/park/victor/scripts/BAMs_txt.py -input_path /n/data1/hms/dbmi/park/victor/other/random/Testfile.txt -output_path /n/data1/hms/dbmi/park/victor/other/tests/ -pipeline Germline -t1 2-00:00:00 -t2 1000 



Actual: 
python3 /n/data1/hms/dbmi/park/victor/scripts/other/BAMs_txt.py -input_path /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/recalibrated_aligned/normal_list -output_path /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/HaplotypeCaller_recalibrated -t1 3-0:00:00 -t2 3000 -r1 1 -r2 2 -input_json /n/data1/hms/dbmi/park/victor/scripts/GATK_Germline_SNPs_Indels/haplotypecaller-gvcf-gatk4.grch37.wgs.inputs.json -queue medium

python3 /n/data1/hms/dbmi/park/victor/scripts/other/BAMs_txt.py -input_path /n/data1/hms/dbmi/park/victor/Doga/Gerburg_WES_2/normal_list -output_path /n/data1/hms/dbmi/park/victor/Doga/Gerburg_WES_2 -t1 4-0:00:00 -t2 3400 -r1 0 -r2 20

python3 /n/data1/hms/dbmi/park/victor/scripts/other/BAMs_txt.py -input_path /n/data1/hms/dbmi/park/ethan/GERBURG/.PreProcessing/bam_list -output_path /n/data1/hms/dbmi/park/ethan/GERBURG/.PreProcessing -t1 5-0:00:00 -t2 3000 -r1 0 -r2 2 

python3 /n/data1/hms/dbmi/park/victor/scripts/other/BAMs_txt.py -input_path /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/bam_list -output_path /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing -t1 3-0:00:00 -queue park -r1 10 -r2 20

python3 /n/data1/hms/dbmi/park/victor/scripts/other/BAMs_txt.py -gatk new -input_path /home/clb36/parkhome/Juan/original/.PreProcessing/bam_list -output_path /home/clb36/parkhome/Juan/original -t1 3-0:00:00 -r1 10 -r2 20 -queue medium

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
    parser.add_argument('-input_path', help='Path to the text file')
    parser.add_argument('-output_path', default='/n/data1/hms/dbmi/park/victor/other/', help='Path to where the final pipeline output will be written to. Default output will be to /n/data1/hms/dbmi/park/victor/other/.')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='0-15:0:0', help='slurm job submission option')
    parser.add_argument('-r1', default=0, help='Lower range bound of indices of BAMs to run')
    parser.add_argument('-r2', default=100000, help='Upper range bound of indices of BAMs to run')
    parser.add_argument('--mem_per_cpu', default='10G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', default=None, help='slurm job submission option')
    parser.add_argument('-gatk', '--gatk_path', default='/home/mk446/BiO/Install/GATK4.1.2.0/gatk', help='path to software execution script')
    parser.add_argument('-reference', '--reference_path', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta', help='path to reference_path file')
    parser.add_argument('-scatter', '--scatter_size', default='50')
    parser.add_argument('-n', '--num_cores', default='1', help='slurm job submission option')

    return parser.parse_args()

def main(): 
    tools_dir = '/n/data1/hms/dbmi/park/victor/scripts/'
    args = vars(parse_args())

    if args['mail_user'] is None:
        print('No email given.')
        sys.exit()

    output_dir, pipeline = arg_clean(args)
    #move_files(args)

    os.system('module load gcc/6.2.0 python/3.6.0 java')

    with open(args['input_path'], 'r') as f:
    	for index, line in enumerate(f):
             if not line.isspace() and index in range(int(args['r1']), int(args['r2'])):
                  os.system('python3 ' + tools_dir + pipeline + ' -in_file ' + line.strip('\n') + ' -out ' + output_dir + ' -t ' + args['runtime'] + ' --mem_per_cpu ' + args['mem_per_cpu']
                        + ' -p ' + args['queue'] + ' --mail_user ' + args['mail_user'] + ' -gatk ' + args['gatk_path'] + ' -reference ' + args['reference_path'] + ' -scatter ' + args['scatter_size'] + 
                        ' -n ' + args['num_cores'])

def arg_clean(args):
    """Cleaning the parsed arguments. Specifically, this function does the following: 
    1. Finds the host filename
    2. Finds the desired output directory
    3. Creates a new output directory within the specified directory
    4. Specifies the appropriate script, depending on the desired pipeline
    """
    if '.txt' in args['input_path']:
        filename = ntpath.basename(args['input_path'])
        filename = filename.split(".")[0]
        #filename = re.findall('/[A-Za-z0-9_]*\.', args['input_path'])[0][1:-1]
    else:
        new_name = args['input_path'] + '.txt'
        filename = ntpath.basename(args['input_path'])
        filename = filename.split(".")[0]
        #filename = re.findall('/[A-Za-z0-9_]*\.', new_name)[0][1:-1]

    if args['output_path'] == '/n/data1/hms/dbmi/park/victor/other/': output_dir = args['output_path']
    if args['output_path'] == '.' or args['output_path'] == './': output_dir = os.getcwd() 
    else: output_dir = args['output_path']

    pipeline = 'GATK_Germline_SNPs_Indels/HaplotypeCaller_4.1.2.0.py'

    return output_dir, pipeline

if __name__ == "__main__":
    main()

