import argparse
import os
import re
import ntpath
import sys
import glob

def parse_args():
    """Uses argparse to enable user to customize script functionality""" 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-tumor', '--input_tumor_path', help='path to input tumor file')
    parser.add_argument('-normal', '--input_normal_path', help='path to normal file')
    parser.add_argument('-csv', default=None, help='csv containing matched tumor/normal pairs')
    parser.add_argument('-out', '--output_directory', default='./', help='directory to which the output directory "/.Mutect2/" will be written to')
    parser.add_argument('-n', '--num_cores', default='1', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='2-00:00:00', help='slurm job submission option')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='15G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', default=None, help='slurm job submission option')
    parser.add_argument('-reference', '--reference_path', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta', help='path to reference_path file')
    # parser.add_argument('-reference', '--reference_path', default='/n/dlsata1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38/Homo_sapiens_assembly38.fasta', help='path to reference_path 
    parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf.gz', help='path to dbsnp file')
    # parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_hg38_common_all_20160601.vcf', help='path to dbsnp file')
    parser.add_argument('-data_type', default='WGS')
    parser.add_argument('-r1', default=1, help='Lower range bound of indices of BAMs to run')
    parser.add_argument('-r2', default=100000, help='Upper range bound of indices of BAMs to run')
    parser.add_argument('-cn', default="1", help='number of cores for Cromwell jobs')
    parser.add_argument('-ct', default="1000", help='cromwell run time; please specify as number of minutes')
    parser.add_argument('-cm', default='7000', help='cromwell cpu memory per core')
    parser.add_argument('-cromwell', '--cromwell_path', default='/n/shared_db/singularity/hmsrc-gatk/cromwell-43.jar', help='path to cromwell.jar file')

    return parser.parse_args()

def arg_clean(args):
    if args['output_directory'] == '.' or args['output_directory'] == './': output_dir = os.getcwd() 
    else: output_dir = args['output_directory']

    return output_dir

def collect_called(normals_path):
    normals = [os.path.realpath(file) for file in glob.glob(os.path.join(normals_path, '*.txt'))]
    return normals

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
    dirname = os.path.dirname(os.path.abspath(__file__))
    tool = os.path.join(dirname, 'MuSE.py')
    args = vars(parse_args())
    if args['mail_user'] is None: 
        print('No email given.')
        sys.exit()
    elif args['csv'] is None:
        print('No csv given.')
        sys.exit()
    output_dir = arg_clean(args)
    

    if int(args['r1']) == 0 and args['mode'].lower() == 'call':
        print("Invalid r1 parameter. r1 must be greater than 0 to account for headers")
        sys.exit()

    os.system('module load gcc/6.2.0 python/3.6.0 java perl')

    tumor_index = get_column(args['csv'], "T")
    normal_index = get_column(args['csv'], "N")
    with open(args['csv'], 'r') as f:
        for index, line in enumerate(f):
            if not line.isspace() and index in range(int(args['r1']), int(args['r2'])):
                current_tumor_sample = get_bam(args['csv'], index, tumor_index)
                current_normal_sample = get_bam(args['csv'], index, normal_index)

                if current_tumor_sample is None or current_normal_sample is None:
                    continue
                else: 
                    tumor_sample = os.path.join(args['input_tumor_path'], get_bam(args['csv'], index, tumor_index))
                    normal_sample = os.path.join(args['input_normal_path'], get_bam(args['csv'], index, normal_index))

                    os.system('python3 ' + tool + ' -tumor ' + tumor_sample + ' -normal ' + normal_sample + ' -out ' + output_dir + ' -t ' + args['runtime'] + ' -n ' + args['num_cores'] + 
                        ' -p ' + args['queue'] + ' --mail_user ' + args['mail_user'] + ' --mem_per_cpu ' + args['mem_per_cpu'] + ' --mail_type ' + args['mail_type'] + ' -reference ' + args['reference_path'] 
                        + ' -data_type ' + args['data_type'] + ' -cn ' + args['cn'] + ' -ct ' + args['ct'] + ' -cm ' + args['cm'] + ' -cromwell ' + args['cromwell_path'] + ' -dbsnp ' + args['dbsnp_path'])


if __name__ == "__main__":
    main()

