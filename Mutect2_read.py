"""This script takes in text files of lists of BAM files and executes the GATK Haplotype and Joint Mutation Callers to produce the desired Output."""

"""
Testing queries:
python3 /n/data1/hms/dbmi/park/victor/scripts/BAMs_txt.py -input_path /n/data1/hms/dbmi/park/victor/other/random/Testfile.txt -output_path /n/data1/hms/dbmi/park/victor/other/tests/ -pipeline Germline
python3 /n/data1/hms/dbmi/park/victor/scripts/BAMs_txt.py -input_path /n/data1/hms/dbmi/park/victor/other/random/Testfile.txt -output_path /n/data1/hms/dbmi/park/victor/other/tests/ -pipeline Germline -t1 2-00:00:00 -t2 1000 



Actual: 
python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2_read.py -tumor /n/data1/hms/dbmi/park/victor/Doga/Gerburg_WES_2 -normal /n/data1/hms/dbmi/park/victor/Doga/Gerburg_WES_2 -csv /n/data1/hms/dbmi/park/victor/Doga/Gerburg_WES_2/Panc_matched_realigned_list.csv -out /n/data1/hms/dbmi/park/victor/Doga/Gerburg_WES_2 -t 4-0:00:00 -r1 40 -p medium 

python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2_read.py -tumor /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK -normal /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK -csv /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK/tumor_bams_cp.csv -out /n/data1/hms/dbmi/park/doga/Gerburg/bam_files_MSK  -t 4-0:00:00 -r1 0 -r2 100 -reference /n/data1/hms/dbmi/park/victor/references/human_g1k_v37.fasta 


python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2_read.py -tumor /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI -normal /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI -csv /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/PancSeq_WES_cohort.csv -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI/Mutect_recalibrated  -t 4-0:00:00 -r2 10

python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2_read.py -tumor /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel -normal /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel -csv /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel/Panc_new_matched_samples.csv -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel -t 3-0:00:00 -reference /n/data1/hms/dbmi/park/victor/references/human_g1k_v37.fasta -r2 2

python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2_read.py -tumor /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel -pon /n/data1/hms/dbmi/park/victor/references/TCGA_1000_PON.b37.vcf -csv /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel/Panc_new_matched_samples.csv -out /n/data1/hms/dbmi/park/DATA/PARP_Panc_DFCI_panel -t 3-0:00:00 -r1 50 -r2 75

python3 /n/data1/hms/dbmi/park/victor/scripts/other/Mutect2_read.py -tumor /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/edited_bams -pon /n/data1/hms/dbmi/park/victor/references/TCGA_1000_PON.hg19.REORDERED.vcf -csv /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/edited_bams/Tumor_Normal.csv -out /n/data1/hms/dbmi/park/doga/FoundationOne/sourceData/Mutect2_4.1.2.0 -t 3-12:00:00 -reference /home/mk446/BiO/Install/GATK-bundle/2.8/hg19/ucsc.hg19.fasta -p medium


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
    parser.add_argument('-pon', '--panel', default='nopath', help='path to panel of norms')
    parser.add_argument('-csv', help='csv containing matched tumor/normal pairs')
    parser.add_argument('-out', '--output_directory', default='./', help='directory to which the output directory "/.Mutect2/" will be written to')
    parser.add_argument('-n', '--num_cores', default='1', help='slurm job submission option')
    parser.add_argument('-t', '--runtime', default='0-12:0:0', help='slurm job submission option')
    parser.add_argument('-p', '--queue', default='park', help='slurm job submission option')
    parser.add_argument('--mem_per_cpu', default='10G', help='slurm job submission option')
    parser.add_argument('--mail_type', default='FAIL', help='slurm job submission option')
    parser.add_argument('--mail_user', default='none', help='slurm job submission option')
    parser.add_argument('-gatk', '--gatk_path', default='/home/mk446/BiO/Install/GATK3.5_160425_g7a7b7cd/GenomeAnalysisTK.jar', help='path to software')
    # parser.add_argument('-gatk', '--gatk_path', default='/n/data1/hms/dbmi/park/lawrence/bin/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar', help='path to software')
    parser.add_argument('-gatk4', '--gatk4_path', default='/n/data1/hms/dbmi/park/alon/software/gatk/gatk-4.0.3.0/gatk', help='path to software')
    parser.add_argument('-reference', '--reference_path', default='/home/mk446/BiO/Install/GATK-bundle/2.8/b37/human_g1k_v37_decoy.fasta', help='path to reference_path file')
    # parser.add_argument('-reference', '--reference_path', default='/n/data1/hms/dbmi/park/SOFTWARE/REFERENCE/hg38/Homo_sapiens_assembly38.fasta', help='path to reference_path 
    parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_b37_common_all_20160601.vcf', help='path to dbsnp file')
    # parser.add_argument('-dbsnp', '--dbsnp_path', default='/home/mk446/BiO/Install/GATK-bundle/dbsnp_147_hg38_common_all_20160601.vcf', help='path to dbsnp file')
    parser.add_argument('-gnomad', '--gnomad_path', default='/n/data1/hms/dbmi/park/victor/software/GATK_bundle/af-only-gnomad.hg19.vcf', help='path to cosmic file' )
    parser.add_argument('-scatter', '--scatter_size', default='50')
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
    if args['output_directory'] == '.' or './': output_dir = os.getcwd() 
    output_dir = os.path.join(args['output_directory'], filename)
    os.makedirs(output_dir, exist_ok = True)
    return output_dir

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
    tools_dir = '/n/data1/hms/dbmi/park/victor/scripts/other/Mutect2.py'
    args = vars(parse_args())
    output_dir = arg_clean(args)

    if int(args['r1']) == 0:
        print("Invalid r1 parameter. r1 must be greater than 0 to account for headers")
        sys.exit()

    os.system('module load gcc/6.2.0 python/3.6.0 java perl')

    #Getting the indices
    tumor_index = get_column(args['csv'], "T")
    normal_index = get_column(args['csv'], "N")

    with open(args['csv'], 'r') as f:
        for index, line in enumerate(f):
            if not line.isspace() and index in range(int(args['r1']), int(args['r2'])):
                current_tumor_sample = get_bam(args['csv'], index, tumor_index)
                current_normal_sample = get_bam(args['csv'], index, normal_index)
                print(current_normal_sample)

                if current_tumor_sample is None or current_normal_sample is None:
                    #continue
                    if current_normal_sample is None and current_tumor_sample is not None and args['panel'] != 'nopath': 
                        tumor_sample = os.path.join(args['input_tumor_path'], get_bam(args['csv'], index, tumor_index))
                        os.system('python3 ' + tools_dir + ' -tumor ' + tumor_sample + ' -pon ' + args['panel'] + ' -out ' + output_dir + ' -t ' + args['runtime'] + ' -n ' + args['num_cores'] +
                            ' -p ' + args['queue'] + ' --mail_user ' + args['mem_per_cpu'] + ' --mail_type ' + args['mail_type'] + ' -reference ' + args['reference_path'] + ' -dbsnp ' + args['dbsnp_path'] + ' -scatter ' + args['scatter_size'] + args['gnomad_path'])
                        #print('python3 ' + tools_dir + ' -tumor ' + tumor_sample + ' -pon ' + args['panel'] + ' -out ' + output_dir + ' -t ' + args['runtime'] + 
                        #    ' -p ' + args['queue'] + ' --mail_user ' + args['mem_per_cpu'] + ' -reference ' + args['reference_path'])                        
                    else:
                        continue
                else: 
                    tumor_sample = os.path.join(args['input_tumor_path'], get_bam(args['csv'], index, tumor_index))
                    normal_sample = os.path.join(args['input_normal_path'], get_bam(args['csv'], index, normal_index))

                    os.system('python3 ' + tools_dir + ' -tumor ' + tumor_sample + ' -normal ' + normal_sample + ' -out ' + output_dir + ' -t ' + args['runtime'] + ' -n ' + args['num_cores'] + 
                        ' -p ' + args['queue'] + ' --mail_user ' + args['mem_per_cpu'] + ' --mail_type ' + args['mail_type'] + ' -reference ' + args['reference_path'] + ' -dbsnp ' + args['dbsnp_path'] + ' -scatter ' + args['scatter_size'] + ' -gnomad ' + args['gnomad_path'])

                    #print(tumor_sample)
                    #print(normal_sample)
                    #print('python3 ' + tools_dir + ' -tumor ' + tumor_sample + ' -normal ' + normal_sample + ' -out ' + output_dir + ' -t ' + args['runtime'] + 
                     #   ' -p ' + args['queue'] + ' --mail_user ' + args['mem_per_cpu'])

if __name__ == "__main__":
    main()

