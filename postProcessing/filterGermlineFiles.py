import argparse
import os
import re
import ntpath
import sys
import csv

def parse_args():
    parser = argparse.ArgumentParser(description='Dataset Creation')

    parser.add_argument('-vcf', type=str, default='ogbl-collab')
    parser.add_argument('-vaf_cut', type=float, default=0)
    parser.add_argument('-ad_cut', type=int, default=0)
    parser.add_argument('-tot_cut', type=int, default=0)
    parser.add_argument('-out', type=str, default=None)
    parser.add_argument('--override', action='store_true')

    return parser.parse_args()

def main(): 
    tools_dir = '/n/data1/hms/dbmi/park/victor/scripts/'
    args = vars(parse_args())
    name = args['vcf'].split(".vcf")[0]
    
    out_file = os.path.join(args['out'], name + ".germline_PASS.vcf")
    filtered_file = os.path.join(args['out'], name + ".germline_filtered.vcf")
    if args['override']: 
        try:
            os.remove(out_file)
        except OSError:
            pass        
        try:
            os.remove(filtered_file)
        except OSError:
            pass
    if not os.path.exists(out_file):
        end = add_vcf_header(args, out_file, True)
        add_vcf_header(args, filtered_file, True)
    
    idx = get_tumor_column(args['vcf'])
    with open(args['vcf'], 'r') as f: 
        for index, line in enumerate(f):
            if not line.isspace() and index >= end: 
                if check_read_levels(line, args, idx): 
                    with open(out_file, 'a+') as o: 
                        o.write(line)
                else:
                    with open(filtered_file, 'a+') as o: 
                        o.write(line) 


def check_read_levels(vcf_line, args, index):
    line_list = vcf_line.rstrip().split('\t')
    #print(len(line_list))
    tumor_info = line_list[index].split(':')[1].split(',')

    if len(tumor_info) != 2: 
        return False
    ref_level = int(tumor_info[0])
    alt_level = int(tumor_info[1])
    if ref_level == 0 and alt_level == 0: 
        return False
    vaf_level = alt_level / (ref_level + alt_level)
    if alt_level < args['ad_cut'] or (ref_level + alt_level) < args['tot_cut'] or vaf_level <= args['vaf_cut']:
        return False
    else:
        return True

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

                for i in range(len(line_list)):
                    if line_list[i] != normal_sample and line_list[i] not in blacklist and gatk_4120: 
                        return i
                    if line_list[i] == 'TUMOR':
                        return i

    return 9

def add_vcf_header(args, file_path, write):
    end = 0
    with open(args['vcf'], 'r') as f_input:
        for index, line in enumerate(f_input):
            if '#' in line:
                if write:
                    with open(file_path, 'a+') as g_out:
                        g_out.write(line)
                end += 1
    return end

if __name__ == "__main__":
    main()
