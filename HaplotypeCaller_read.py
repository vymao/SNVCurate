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
    parser.add_argument('-pipeline', default='Germline', help='The corresponding pipeline for the method calls (ie. Somatic or Germline)')
    parser.add_argument('-queue', default='park', help='The queue used for running jobs')
    parser.add_argument('-t1', default='5-00:00:00', help='HaplotypeCaller runtime')
    parser.add_argument('-t2', default='2600', help='Cromwell task runtime')
    parser.add_argument('-r1', default=0, help='Lower range bound of indices of BAMs to run')
    parser.add_argument('-r2', default=100000, help='Upper range bound of indices of BAMs to run')
    parser.add_argument('-m', default='victor_mao@hms.harvard.edu', help='Email for notifications')
    parser.add_argument('-gatk', default='old', help='GATK version: old or new')
    parser.add_argument('-input_json', default='/n/data1/hms/dbmi/park/victor/scripts/GATK_Germline_SNPs_Indels/haplotypecaller-gvcf-gatk4.hg37.wgs.inputs.json', help='Reference genome-dependent. Path to gatk4-haplotypecaller file')

    return parser.parse_args()

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
    output_dir = os.path.join(args['output_path'], filename + '.HaplotypeCaller')
    os.makedirs(output_dir, exist_ok = True)

    if args['pipeline'].lower() == 'somatic':
        pipeline = 'GATK_Somatic_SNPs_Indels/Mutect2.py'
    else:
        if args['gatk'] == 'new':
            pipeline = 'GATK_Germline_SNPs_Indels/HaplotypeCaller_4.1.2.0.py'
        elif args['queue'] == 'priopark' or args['queue'] == 'park':
            pipeline = 'GATK_Germline_SNPs_Indels/HaplotypeCaller_original.py'
        else:
            pipeline = 'GATK_Germline_SNPs_Indels/HaplotypeCaller.py'

    return output_dir, pipeline

def move_files(args):
    """Moves files to a semi-temporary directory under /n/scratch2/vym1 on the O2 cluster to group all files in a common directory"""
    home = '/n/scratch2/vym1'
    filename = re.findall('/[A-Za-z0-9_]*\.', args['input_path'])[0][1:-1]
    newhome = os.path.join(home, filename)
    temp_out = os.makedirs(newhome, exist_ok = True)

    with open(args['input_path'], 'r') as f:
        for line in f:
            if not line.isspace():
                parent_directory = re.split('^(.+)/([^/]+)$', line)[1]
                filename = re.findall('/[A-Za-z0-9_]*\.', line)[0][1:-1]
                for file in os.listdir(parent_directory):
                    if filename in file:
                        cmd = ['sbatch /n/data1/hms/dbmi/park/victor/other/random/Move_file.sh', os.path.join(parent_directory, file.strip('\n')), newhome]
                        subprocess.Popen(cmd).wait()
                        #os.system('sbatch /n/data1/hms/dbmi/park/victor/other/random/Move_file.sh ' + os.path.join(parent_directory, file.strip('\n')) + ' ' + newhome)

def main(): 
    tools_dir = '/n/data1/hms/dbmi/park/victor/scripts/'
    args = vars(parse_args())
    output_dir, pipeline = arg_clean(args)
    #move_files(args)

    os.system('module load gcc/6.2.0 python/3.6.0 java')
    #temp_out = '/n/scratch2/vym1'  

    if args['pipeline'].lower() == 'somatic':
        with open(args['input_path'], 'r') as f:
            for index, line in enumerate(f):
                if not line.isspace() and index in range(int(args['r1']), int(args['r2'])):
                    os.system('python3 ' + tools_dir + pipeline + ' -in_file ' + line.strip('\n') + ' -out ' + output_dir + ' -t ' + args['t1'] + ' -t2 ' + args['t2']) 
    else: 
        with open(args['input_path'], 'r') as f:
            for index, line in enumerate(f):
                if not line.isspace() and index in range(int(args['r1']), int(args['r2'])):
                    #print('python3 ' + tools_dir + pipeline + ' -in_file ' + line.strip('\n') + ' -out ' + output_dir + ' -t ' + args['t1'] + ' -t2 ' + args['t2'] + ' -p ' + args['queue'] + ' --mail_user ' + args['m'])
                    os.system('python3 ' + tools_dir + pipeline + ' -in_file ' + line.strip('\n') + ' -out ' + output_dir + ' -t ' + args['t1'] + ' -t2 ' + args['t2'] + ' -p ' + args['queue'] + ' --mail_user ' + args['m'])
if __name__ == "__main__":
    main()

