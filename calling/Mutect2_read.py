import argparse
import os
import re
import ntpath
import sys

def parse_args():
    """Uses argparse to enable user to customize script functionality""" 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-tumor', '--input_tumor_path', help='path to input tumor file')
    parser.add_argument('-normal', '--input_normal_path', help='path to normal file')
    parser.add_argument('-pon', '--panel', default='nopath', help='path to panel of norms')
    parser.add_argument('-csv', help='csv containing matched tumor/normal pairs')
    parser.add_argument('-out', '--output_directory', default='./', help='directory to which the output directory "/.Mutect2/" will be written to')
    parser.add_argument('-n', '--num_cores', default='1', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='0-12:0:0', help='slurm job submission option')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='10G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', default=None, help='slurm job submission option')
    parser.add_argument('-reference', '--reference_path', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta', help='path to reference_path file')
    # parser.add_argument('-reference', '--reference_path', default='/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38/Homo_sapiens_assembly38.fasta', help='path to reference_path 
    parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf', help='path to dbsnp file')
    # parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_hg38_common_all_20160601.vcf', help='path to dbsnp file')
    parser.add_argument('-gnomad', '--gnomad_path', default='/n/data1/hms/dbmi/park/victor/software/GATK_bundle/af-only-gnomad.raw.sites.b37.vcf.gz', help='path to gnomad file' )
    parser.add_argument('-scatter', '--scatter_size', default='50')
    parser.add_argument('-r1', default=1, help='Lower range bound of indices of BAMs to run')
    parser.add_argument('-r2', default=100000, help='Upper range bound of indices of BAMs to run')
    parser.add_argument('-cn', default="1", help='number of cores for Cromwell jobs')
    parser.add_argument('-ct', default="1000", help='cromwell run time; please specify as number of minutes')
    parser.add_argument('-cm', default='7000', help='cromwell cpu memory per core')
    parser.add_argument('-cromwell', default='/n/shared_db/singularity/hmsrc-gatk/cromwell-43.jar', help='path to cromwell.jar file')
    parser.add_argument('-interval_list', default='/n/data1/hms/dbmi/park/victor/software/MuTecT2_b37_scattered_intervals.txt')
    parser.add_argument('-gatk', '--gatk_path', default='/home/mk446/BiO/Install/GATK4.1.2.0/gatk', help='path to software')
    parser.add_argument('-gatk4', '--gatk4_path', default='/n/data1/hms/dbmi/park/alon/software/gatk/gatk-4.0.3.0/gatk', help='path to software')
    parser.add_argument('-parallel', default='True')
    return parser.parse_args()

def arg_clean(args):
    if args['output_directory'] == '.' or args['output_directory'] == './': output_dir = os.getcwd()
    else: output_dir = args['output_directory'] 

    return output_dir

def generate_regions_files(args):
    regions_out_directory = os.path.join(args['output_directory'], '.Mutect2/.regions/')
    os.makedirs(regions_out_directory, exist_ok=True)
    os.system(args['gatk4_path'] + ' SplitIntervals' + '\\' + '\n' + \
     '\t' + '-R ' + args['reference_path'] + ' \\' + '\n' + \
     '\t' + '-scatter ' + args['scatter_size'] + ' \\' + '\n' + \
     '\t' + '-O ' + regions_out_directory + ' \\')

def return_region_files(args):
    regions_out_directory = os.path.join(args['output_directory'], '.Mutect2/.regions/')
    region_files = [os.path.join(regions_out_directory, file) for file in os.listdir(regions_out_directory) if "scattered.interval_list" in file]
    return region_files

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
    tool = os.path.join(dirname, 'Mutect2.py')
    args = vars(parse_args())
    output_dir = arg_clean(args)

    if args['mail_user'] is None:
        print('No email given.')
        sys.exit()

    if int(args['r1']) == 0:
        print("Invalid r1 parameter. r1 must be greater than 0 to account for headers")
        sys.exit()

    os.system('module load gcc/6.2.0 python/3.6.0 java bcftools')
    
    reference_name = os.path.basename(args['reference_path']).split('.')[0]
    if (args['scatter_size'] != '50' or args['reference_path'] != '/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta') and args['parallel'].lower() == "true":
        regions_out_directory = os.path.join(args['output_directory'], '.Mutect2/.regions/')
        intervals_list_file = os.path.join(regions_out_directory, reference_name + '.scattered_intervals.list')
        if os.path.isfile(intervals_list_file):
            os.remove(intervals_list_file) 
        generate_regions_files(args)
        region_files = return_region_files(args)
        with open(intervals_list_file, 'a') as intervals:
            for file in region_files:
                intervals.write(file + '\n')
        args['interval_list'] = intervals_list_file

    tumor_index = get_column(args['csv'], "T")
    normal_index = get_column(args['csv'], "N")

    with open(args['csv'], 'r') as f:
        for index, line in enumerate(f):
            if not line.isspace() and index in range(int(args['r1']), int(args['r2'])):
                current_tumor_sample = get_bam(args['csv'], index, tumor_index)
                current_normal_sample = get_bam(args['csv'], index, normal_index)

                if current_tumor_sample is None or current_normal_sample is None:
                    if current_normal_sample is None and current_tumor_sample is not None and args['panel'] != 'nopath': 
                        tumor_sample = os.path.join(args['input_tumor_path'], get_bam(args['csv'], index, tumor_index))
                        os.system('python3 ' + tool + ' -tumor ' + tumor_sample + ' -pon ' + args['panel'] + ' -out ' + output_dir + ' -t ' + args['runtime'] + ' -n ' + args['num_cores'] +
                            ' -p ' + args['queue'] + ' --mail_user ' + args['mail_user'] + ' --mem_per_cpu ' + args['mem_per_cpu'] + ' --mail_type ' + args['mail_type'] + ' -reference ' + args['reference_path'] 
                            + ' -dbsnp ' + args['dbsnp_path'] + ' -scatter ' + args['scatter_size'] + ' -gnomad ' + args['gnomad_path'] + ' -cn ' + args['cn'] + ' -ct ' + args['ct'] + ' -cm ' + args['cm'] 
                            + ' -cromwell ' + args['cromwell'] + ' -parallel ' + args['parallel'])                    
                    else:
                        continue
                else: 
                    tumor_sample = os.path.join(args['input_tumor_path'], get_bam(args['csv'], index, tumor_index))
                    normal_sample = os.path.join(args['input_normal_path'], get_bam(args['csv'], index, normal_index))

                    os.system('python3 ' + tool + ' -tumor ' + tumor_sample + ' -normal ' + normal_sample + ' -out ' + output_dir + ' -t ' + args['runtime'] + ' -n ' + args['num_cores'] + 
                        ' -p ' + args['queue'] + ' --mail_user ' + args['mail_user'] + ' --mem_per_cpu ' + args['mem_per_cpu'] + ' --mail_type ' + args['mail_type'] + ' -reference ' + args['reference_path'] 
                        + ' -dbsnp ' + args['dbsnp_path'] + ' -scatter ' + args['scatter_size'] + ' -gnomad ' + args['gnomad_path'] + ' -cn ' + args['cn'] + ' -ct ' + args['ct'] + ' -cm ' + args['cm'] 
                        + ' -cromwell ' + args['cromwell'] + ' -parallel ' + args['parallel'] + ' -interval_list ' + args['interval_list'])

if __name__ == "__main__":
    main()

