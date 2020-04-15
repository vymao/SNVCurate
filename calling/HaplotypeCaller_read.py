import argparse
import os
import re
import ntpath
import sys

def parse_args():
    """Uses argparse to enable user to customize script functionality""" 
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
    parser.add_argument('-n', '--num_cores', default='1', help='slurm job submission option')
    parser.add_argument('-cn', default="1", help='number of cores for Cromwell jobs')
    parser.add_argument('-ct', default="3000", help='cromwell run time; please specify as number of minutes')
    parser.add_argument('-cm', default='5000', help='cromwell cpu memory per core')
    parser.add_argument('-cromwell', '--cromwell_path', default='/n/shared_db/singularity/hmsrc-gatk/cromwell-43.jar', help='path to cromwell.jar file')
    parser.add_argument('-reference_name', default='b37', help='hg19, b37, etc.')

    return parser.parse_args()

def main(): 
    dirname = os.path.dirname(os.path.abspath(__file__))
    tool = os.path.join(dirname, 'HaplotypeCaller.py')
    args = vars(parse_args())

    if args['mail_user'] is None:
        print('No email given.')
        sys.exit()

    output_dir = arg_clean(args)

    os.system('module load gcc/6.2.0 python/3.6.0 java')

    with open(args['input_path'], 'r') as f:
    	for index, line in enumerate(f):
             if not line.isspace() and index in range(int(args['r1']), int(args['r2'])):
                  os.system('python3 ' + tool + ' -in_file ' + line.strip('\n') + ' -out ' + output_dir + ' -t ' + args['runtime'] + ' --mem_per_cpu ' + args['mem_per_cpu']
                        + ' -p ' + args['queue'] + ' --mail_user ' + args['mail_user'] + ' -gatk ' + args['gatk_path'] + ' -reference ' + args['reference_path'] + ' -n ' + args['num_cores']
                        + ' -cn ' + args['cn'] + ' -ct ' + args['ct'] + ' -cm ' + args['cm'] + ' -cromwell ' + args['cromwell'] + ' -reference_name ' + args['reference_name'])

def arg_clean(args):
    if args['output_path'] == '/n/data1/hms/dbmi/park/victor/other/': output_dir = args['output_path']
    if args['output_path'] == '.' or args['output_path'] == './': output_dir = os.getcwd() 
    else: output_dir = args['output_path']

    return output_dir

if __name__ == "__main__":
    main()

