import argparse
import os
import re
import ntpath
import sys
import csv

def parse_args():
    """Uses argparse to enable user to customize script functionality""" 
    #BE SURE TO ADD MULTIPLE TEXT FILE ARGUMENTS Later
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-input_path', help='Path to the input file')
    parser.add_argument('-vcf_path', help='Path to vcf file')
    parser.add_argument('-output_path', default='/n/data1/hms/dbmi/park/victor/other/', help='Path to where the final pipeline output will be written to. Default output will be to /n/data1/hms/dbmi/park/victor/other/.')
    parser.add_argument('-file_type', default='anno', help='Type of file for cut: vcf or annotated csv (anno)')
    parser.add_argument('-queue', default='park', help='The queue used for running jobs')
    parser.add_argument('-t1', default='5-00:00:00', help='HaplotypeCaller runtime')
    parser.add_argument('-t2', default='2600', help='Cromwell task runtime')
    parser.add_argument('-r1', default=0, help='Lower range bound of indices of BAMs to run')
    parser.add_argument('-r2', default=100000, help='Upper range bound of indices of BAMs to run')
    parser.add_argument('-m', default='victor_mao@hms.harvard.edu', help='Email for notifications')
    parser.add_argument('-cut', '--threshhold', default=0.01, help='Threshhold cutoff value for annotations')
    parser.add_argument('-mode', default='annotate', help='Mode of operation: annotate or combine')
    parser.add_argument('-hap_path', help='Path to HaplotypeCaller annotated output')
    parser.add_argument('-hap', default=None, help='Path to matched normal file')
    parser.add_argument('-mut', help='Path to Mutect2 subdirectories')
    parser.add_argument('-normals', help='Path to matched normals csv')
    parser.add_argument('-pon', help='Panel of normals for samples without normals')
    parser.add_argument('-alt_cut', default=8, help='Alternate read level cutoff')
    parser.add_argument('-tot_cut', default=20, help='Alternate + reference total read level cutoff')
    parser.add_argument('-vaf_cut', default=0.05, help='VAF level cutoff')
    parser.add_argument('-panel', default=False, help='Used Panel of Normals?')

    return parser.parse_args()

def main(): 
    tools_dir = '/n/data1/hms/dbmi/park/victor/scripts/'
    args = vars(parse_args())
    output_dir, base_name = arg_clean(args)


    os.system('module load gcc/6.2.0 python/3.6.0 java') 

    levels = str(args['alt_cut']) + '_' + str(args['tot_cut']) + '_' + str(args['vaf_cut'])

    if args['mode'].lower() == 'annotate':
        if args['file_type'].lower() == 'anno':
            germline_file_path = os.path.join(output_dir, base_name + '.ANNO.germline_variants_filtered.' + levels + '.vcf')
            somatic_file_path = os.path.join(output_dir, base_name + '.ANNO.somatic_variants_filtered.' + levels + '.vcf')
            bad_read_path = os.path.join(output_dir, base_name + '.ANNO.bad_somatic_quality.' + levels + '.vcf')
        else: 
            germline_file_path = os.path.join(output_dir, base_name + '.M2_RISK.germline_variants_filtered.vcf')
            somatic_file_path = os.path.join(output_dir, base_name + '.M2_RISK.somatic_variants_filtered.vcf')

        if not os.path.exists(germline_file_path):
            end = add_vcf_header(args, germline_file_path, True)
        if not os.path.exists(somatic_file_path) and args['file_type'].lower() == 'anno':
            end = add_vcf_header(args, somatic_file_path, True)
        if args['file_type'].lower() == 'anno':
            if not os.path.exists(bad_read_path): 
                end = add_vcf_header(args, bad_read_path, True)

        else: 
            end = add_vcf_header(args, germline_file_path, False)


        if args['file_type'].lower() == 'vcf':
            with open(args['input_path'], 'r') as f:
                for index, line in enumerate(f):
                    if not line.isspace() and index in range(int(args['r1']) + end, int(args['r2']) + end):
                        if not args['panel']:
                            vcf_examine(line, args, germline_file_path, somatic_file_path, 'germline')
                        else: 
                            vcf_examine(line, args, germline_file_path, somatic_file_path, 'panel_of_normals')
        else: 
            if args['hap'] is not None: 
                args['hap_path'] = args['hap']
            else: 
                normal = get_normal(args)
                if normal is None: 
                    args['hap_path'] = args['pon']
                else: 
                    normal_vcf = normal.split('.bam')[0] + '.vcf'
                    args['hap_path'] = os.path.join(args['hap_path'], normal_vcf)

            tumor_index = get_tumor_column(args['vcf_path'])
            with open(args['input_path'], 'r') as f:
                for index, line in enumerate(f):
                    if not line.isspace() and index in range(int(args['r1']) + 1, int(args['r2']) + 1):
                        anno_examine(line, args, germline_file_path, somatic_file_path, bad_read_path, tumor_index)
    else: 
        print('Invalid mode.')

def arg_clean(args):
    """Cleaning the parsed arguments. Specifically, this function does the following: 
    1. Finds the host filename
    2. Finds the desired output directory
    3. Creates a new output directory within the specified directory
    4. Specifies the appropriate script, depending on the desired pipeline
    """
    if args['output_path'] == '/n/data1/hms/dbmi/park/victor/other/': output_dir = args['output_path']
    if args['output_path'] == '.' or args['output_path'] =='./': output_dir = os.getcwd() 
    else: output_dir = args['output_path']

    if args['mode'].lower() == 'annotate':
        filename = ntpath.basename(args['input_path'])
        if args['file_type'].lower() == 'vcf':
            if '.vcf' in args['input_path']: 
                base_name = filename.split(".vcf")[0]
            else:
                base_name = filename
        elif args['file_type'].lower() == 'anno':
            if '.csv' in args['input_path']:
                base_name = filename.split('.csv')[0]
            else:
                base_name = filename
        else:
            print('Invalid file type choice. See options for more details.')
            sys.exit()
    else: 
        base_name = None

    args['alt_cut'] = int(args['alt_cut'])
    args['tot_cut'] = int(args['tot_cut'])
    args['vaf_cut'] = float(args['vaf_cut'])

    return output_dir, base_name

def get_normal(args): 
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

    def get_bam(csv, tumor, column, tumor_column):
        with open(csv, 'r') as file:
            for idx, line in enumerate(file):
                line_list = line.strip().split(',')
                tumor_sample = line_list[tumor_column].split('.')[0]
                if tumor_sample in tumor: 
                    result = line_list[column]
                    if result != 'N/A':
                        return result


    tumor_index = get_column(args['normals'], "T")
    normal_index = get_column(args['normals'], "N")
    current_normal_sample = get_bam(args['normals'], ntpath.basename(args['input_path']), normal_index, tumor_index)
    return current_normal_sample

def vcf_examine(line, args, germline_file_path, somatic_file_path, examine_phrase):
    if examine_phrase in line: 
        with open(germline_file_path, 'a') as f:
            f.write(line)


def anno_examine(line, args, germline_file_path, somatic_file_path, bad_read_path, tumor_index): 
    line_list = line.strip().split('\t')
    germline = False
    #print(line_list[4])
    #print(line_list[27])
    #print("")
    columns_examine = get_annotated_columns(args['input_path'])
    if line_list[1] == "41256830": 
        print('yes')
    
    for item in line_list: 
        if line_list.index(item) not in columns_examine: 
            continue
        if line_list.index(item) > columns_examine[0]:
            try: 
                number = float(item)
                if number > float(args['threshhold']):
                    vcf_line = find_in_vcf(args, line, 'vcf_path')
                    with open(germline_file_path, 'a') as f:
                        if not check_in_file(vcf_line, germline_file_path):
                            f.write(vcf_line)
                            germline = True
                            break
            except:
                continue
            if germline:
                return

    if not germline: 
        vcf_line = find_in_vcf(args, line, 'vcf_path')
        with open(somatic_file_path, 'a') as f:
            
            if check_germline_risk(vcf_line) and find_in_vcf(args, line, 'hap_path') is None and vcf_line is not None:
                if check_read_levels(vcf_line, args, tumor_index): 
                    f.write(vcf_line)
                    return 
                else: 
                    with open(bad_read_path, 'a') as bad: 
                        bad.write(vcf_line)
                        return
    

def check_germline_risk(vcf_line):
    if 'germline_risk' in vcf_line: 
        return False
    else: 
        return True

def check_read_levels(vcf_line, args, index):
    line_list = vcf_line.rstrip().split('\t')
    print(index)
    print(line_list)
    tumor_info = line_list[index].split(':')[1].split(',')

    if len(tumor_info) != 2: 
        return False
    ref_level = int(tumor_info[0])
    alt_level = int(tumor_info[1])
    if ref_level == 0 and alt_level == 0: 
        return False
    vaf_level = alt_level / (ref_level + alt_level)
    if alt_level < args['alt_cut'] or (ref_level + alt_level) < args['tot_cut'] or vaf_level <= args['vaf_cut']:
        return False
    else:
        return True

def check_in_file(search, file): 
    with open(file, 'r') as f:
        for line in f:
            if line == search:
                return True
            else:
                return False

def find_in_vcf(args, line, file):
    line_list = line.rstrip().split('\t')

    chrom = line_list[0]
    pos = line_list[1]

    with open(args[file], 'r') as vcf: 
        for index, vcf_line in enumerate(vcf):
            vcf_line_list = vcf_line.rstrip().split('\t')
            if '#' in line: 
                continue
            elif line.isspace() or chrom != vcf_line_list[0]:
                continue
            else:
                #print(vcf_line_list)
                if int(pos) == int(vcf_line_list[1]) or int(pos) == int(vcf_line_list[1]) - 1 or int(pos) == int(vcf_line_list[1]) + 1:
                    return vcf_line

def add_vcf_header(args, file_path, write):
    end = 0
    with open(args['vcf_path'], 'r') as f_input:
        for index, line in enumerate(f_input):
            if '#' in line:
                if write:
                    with open(file_path, 'a+') as g_out:
                        g_out.write(line)
                end += 1
    return end

def get_tumor_column(file_path):
    blacklist = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    normal_sample = ""
    normal_column = 0
    gatk_4120 = False
    with open(file_path, 'r') as vcf: 
        for index, line in enumerate(vcf):
            if "##normal_sample" in line: 
                normal_sample = line.rstrip().split('=')[1]
                gatk_4120 = True
            if "#CHROM" in line: 
                line_list = line.rstrip().split('\t')
                print(line_list)
                for i in range(len(line_list)):
                    if line_list[i] != normal_sample and line_list[i] not in blacklist and gatk_4120: 
                        return i
                    if line_list[i] == 'TUMOR':
                        return i

    return 9


def get_annotated_columns(file): 
    list_of_databases = ['ExAC', 'gnomAD', '1000g']
    header = ""
    column_indices = []

    with open(file, 'r') as main: 
        for index, line in enumerate(main): 
            if index == 0: 
                header = line
                break

    line_list = line.rstrip().split('\t')
    for item in line_list: 
        for database in list_of_databases: 
            if database in item: 
                column_indices.append(line_list.index(item))

    return column_indices

if __name__ == "__main__":
    main()

