"""This script takes in text files of lists of BAM files and executes the GATK Haplotype and Joint Mutation Callers to produce the desired Output."""

"""
Actual: 
python3 /n/data1/hms/dbmi/park/victor/scripts/other/FMG_Label.py -source Mutect -out /n/data1/hms/dbmi/park/ethan/GERBURG/.PreProcessing/Mutect2_Recalibrated/.Mutect2/Sample_DS-bkm-086-T_RECAL -in /n/data1/hms/dbmi/park/ethan/GERBURG/.PreProcessing/Mutect2_Recalibrated/.Mutect2/Sample_DS-bkm-086-T_RECAL/Sample_DS-bkm-086-T_RECAL.PASS.ANNO.germline_variants_filtered.hg19_multianno.txt.csv.txt


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
    parser.add_argument('-in', '--input_file', help='Path to the input file')
    parser.add_argument('-out', help='Output directory')
    parser.add_argument('-source', help='Source of input file')

    return parser.parse_args()

def get_file(path, end, contains): 
    file_list = []
    for root, direct, file in os.walk(path):
        for file_name in file:
            rel_dir = os.path.relpath(root, path)
            rel_file = os.path.join(rel_dir, file_name)
            if file_name.lower().endswith(end) and contains in file_name.lower():
                file_list.append(os.path.join(path, rel_file,))

    return file_list

def get_last_item(args): 
    with open(args['input_file'], 'r') as file: 
        for index, line in enumerate(file): 
            if '.txt' in args['input_file']: 
                #print('yes')
                line_list = line.strip().split('\t')
            else: 
                line_list = line.rstrip().split(',')            
            if index == 0: 
                header_line_len = len(line_list)
                body_line_len = header_line_len
                #print(header_line_len)
            elif index == 1:   
                body_line_len = len(line_list)
                #print(body_line_len)
            else: 
                #print('here')
                #print(body_line_len - header_line_len)
                return body_line_len - header_line_len, body_line_len
    return body_line_len - header_line_len, body_line_len

def main(): 
    tools_dir = '/n/data1/hms/dbmi/park/victor/scripts/other/Filter_Mutect_Germlines.py'
    args = vars(parse_args())

    os.system('module load gcc/6.2.0 python/3.6.0 java')

    if args['source'].lower() == 'mutect':
        end = 'Germline_Risk'
    elif args['source'].lower() == 'haplotypecaller':
        end = 'Haplotypecaller'
    elif args['source'].lower() == 'filtering':
        end = 'Common variant filtering'
    elif args['source'].lower() == 'haplotypecaller_panel':
        end = 'HaplotypeCaller_PANEL'

    #files = get_file(args['input_file'], 'combined.vcf', 'combined') 
    base = os.path.basename(args['input_file'])
    out = os.path.join(args['out'], base)

    output_name = out + '.LABELED'
    
    line_diff, body_line_len = get_last_item(args)
    total_len = body_line_len

    
    
    with open(args['input_file'], 'r') as file: 
        for index, line in enumerate(file): 
            if '.txt' in args['input_file']: 
                line_list = line.strip().split('\t')
                line_len = len(line_list)
            else: 
                line_list = line.rstrip().split(',')
                line_len = len(line_list)
            with open(output_name, 'a+') as new:
                if index == 0:
                    '''
                    for i in range(line_diff): 
                        line_list.append('.')
                    '''
                    for i in range(total_len - line_len): 
                        line_list.append('.')
                    line_list.append('Source')
                    if '.txt' in args['input_file']:
                        str_data = '\t'.join(str(i) for i in line_list)
                    else: 
                        str_data = ','.join(str(i) for i in line_list)
                    str_data += '\n'
                    new.write(str_data)
                else: 
                    for i in range(total_len - line_len): 
                        line_list.append('.')
                    line_list.append(end)
                    if '.txt' in args['input_file']:
                        str_data = '\t'.join(str(i) for i in line_list)
                    else: 
                        str_data = ','.join(str(i) for i in line_list)  
                    str_data += '\n'  
                    new.write(str_data)   
               

if __name__ == "__main__":
    main()

