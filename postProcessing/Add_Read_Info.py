import argparse
import os
import re
import ntpath
import sys

def parse_args():
    """Uses argparse to enable user to customize script functionality""" 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-in_file', help='Path to the input file')
    parser.add_argument('-vcf_path', help='Path to vcf file')
    parser.add_argument('-hap', default=False, help='HaplotypCaller output?')
    parser.add_argument('-vaf_filter', type=float, default=0)
    parser.add_argument('-ad_filter', type=int, default=0)
    parser.add_argument('-tot_filter', type=int,  default=0)

    return parser.parse_args()

def find_in_vcf(args, line):
    if '.txt' in args['in_file']: 
        line_list = line.strip().split('\t')
    else: 
        lie_list = line.rstrip().split(',')

    chrom = line_list[0]
    pos = line_list[1]

    with open(args['vcf_path'], 'r') as vcf: 
        for index, vcf_line in enumerate(vcf):
            vcf_line_list = vcf_line.rstrip().split('\t')
            if '#' in line: 
                continue
            elif line.isspace() or chrom != vcf_line_list[0]:
                continue
            else:
                if int(pos) == int(vcf_line_list[1]) or int(pos) == int(vcf_line_list[1]) - 1 or int(pos) == int(vcf_line_list[1]) + 1:
                    return vcf_line

def get_read_levels(args, vcf_line, index):
    if vcf_line is None: 
    	return
    line_list = vcf_line.rstrip().split('\t')

    tumor_info = line_list[index].split(':')[1].split(',')

    if len(tumor_info) != 2: 
        return 
    ref_level = int(tumor_info[0])
    alt_level = int(tumor_info[1])
    if ref_level == 0 and alt_level == 0: 
        return [0, 0, 0]
    vaf_level = alt_level / (ref_level + alt_level)
    return [ref_level, alt_level, vaf_level]


def get_last_item(args): 
    with open(args['in_file'], 'r') as file: 
        for index, line in enumerate(file): 
            if '.txt' in args['in_file']: 
                line_list = line.strip().split('\t')
            else: 
                line_list = line.rstrip().split(',')            
            if index == 0: 
                header_line_len = len(line_list)
                body_line_len = header_line_len
            elif index == 1:   
                body_line_len = len(line_list)
            else: 
                return body_line_len - header_line_len, body_line_len
    return (body_line_len - header_line_len, body_line_len)

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

def main(): 
    args = vars(parse_args())
    os.system('module load gcc/6.2.0 python/3.6.0 java')

    new_file = args['in_file'] + '.levels'
    line_diff, body_line_len = get_last_item(args)
    total_len = body_line_len
    tumor_index = get_tumor_column(args['vcf_path'])


    with open(args['in_file'], 'r') as in_file:
        for index, line in enumerate(in_file):
            with open(new_file, 'a+') as new:           
                if index == 0:
                    if '.txt' in args['in_file']: 
                        header = line.strip().split('\t')
                        header_len = len(header)
                    else: 
                        header = line.rstrip().split(',')
                        header_len = len(header)

                    for i in range(total_len - header_len): 
                        header.append('.')

                    header.append('RD')
                    header.append('AD')
                    header.append('vaf_level')
                    if '.txt' in args['in_file']:
                        str_data = '\t'.join(str(i) for i in header)
                    else: 
                        str_data = ','.join(str(i) for i in header)
                    str_data += '\n'
                    new.write(str_data)

                elif '#' not in line: 
                    if '.txt' in args['in_file']: 
                        vcf_line_list = line.strip().split('\t')
                        line_len = len(vcf_line_list)
                    else: 
                        vcf_line_list = line.rstrip().split(',')
                        line_len = len(vcf_line_list)
                    if not args['hap']:
                        for i in range(total_len - line_len): 
                            vcf_line_list.append('.')
    
                        new_line = find_in_vcf(args, line)
                        read_levels = get_read_levels(args, new_line, tumor_index)
                        if read_levels is not None and read_levels[2] >= args['vaf_filter'] and read_levels[1] >= args['ad_filter'] and read_levels[0] + read_levels[1] >= args['tot_filter']: 
                            vcf_line_list.append(read_levels[0])
                            vcf_line_list.append(read_levels[1])
                            vcf_line_list.append(read_levels[2])
                            if '.txt' in args['in_file']:
                                str_data = '\t'.join(str(i) for i in vcf_line_list)
                            else: 
                                str_data = ','.join(str(i) for i in vcf_line_list)
                            str_data += '\n'
                            new.write(str_data)
                    else: 
                        for i in range(total_len - line_len): 
                            vcf_line_list.append('.')
                        vcf_line_list.append('.')
                        vcf_line_list.append('.')
                        vcf_line_list.append('.')
                        if '.txt' in args['in_file']:
                            str_data = '\t'.join(str(i) for i in vcf_line_list)
                        else: 
                            str_data = ','.join(str(i) for i in vcf_line_list)
                        str_data += '\n'
                        new.write(str_data)


if __name__ == "__main__":
    main()

