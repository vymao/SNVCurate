import pandas as pd
import numpy as np
import glob
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser() 
    parser.add_argument('-input_dir')
    parser.add_argument('-csv')
    parser.add_argument('-cut')
    parser.add_argument('-out')

    return parser.parse_args()

def main():
    args = vars(parse_args())
    path_to_coverage = args['input_dir']
    coverage_exclusive_path = os.path.join(path_to_coverage, 'coverage_exclusive') 

    new_files = []
    with open(args['csv'], 'r') as file:
        tumor_index = 0
        for index, line in enumerate(file): 
            line_split = line.rstrip().split(',')
            if index == 0: 
                if line_split[1] == 'Tumor': tumor_index = 1
            else: 
                tumor_sample = line_split[tumor_index]
                sample = tumor_sample.split('.')[0]
                coverage_sample = sample + '.coverage'
                new_files.append(os.path.join(path_to_coverage, coverage_sample))


    for file in new_files: 
        sample = os.path.basename(file).split('.')[0]

        df = pd.read_csv(file, sep='\t', names = ['CHROM', 'START', 'END', 'COVERAGE'])
        chromosomes = df['CHROM'].unique().tolist()
        chromosomes = [str(x) for x in chromosomes]
        chromosomes = [x for x in chromosomes if 'GL' not in x]
        chromosomes = [x for x in chromosomes if 'MT' not in x]
        chromosomes = [x for x in chromosomes if 'NC' not in x]
        df['CHROM'] = df['CHROM'].astype(str)

        coverage_chr = {}
        for chrm in chromosomes: 
            coverage_chr[chrm] = df[df['CHROM'] == chrm]

        test = pd.DataFrame([[]])
        for chrom in coverage_chr.keys():
            filtered = coverage_chr[chrom]
            filtered = filtered[filtered['COVERAGE'] >= int(args['cut'])].drop(columns = ['COVERAGE'])
            test = pd.concat([test, filtered])

        cov_cut_bed = test
        df2 = pd.DataFrame([[]])
        current_chrom = cov_cut_bed.iloc[0][0]
        current_start = cov_cut_bed.iloc[0][1]
        current_end = cov_cut_bed.iloc[0][2]
        for i in range(1, len(cov_cut_bed)):
            if current_end == cov_cut_bed.iloc[i][1]:
                current_end = cov_cut_bed.iloc[i][2]
            else: 
                grouped_bin = pd.DataFrame([[current_chrom, current_start, current_end]])
                df2 = pd.concat([df2, grouped_bin], ignore_index = True)
                current_chrom = cov_cut_bed.iloc[i][0]
                current_start = cov_cut_bed.iloc[i][1]
                current_end = cov_cut_bed.iloc[i][2]
                
        df2 = df2.dropna()
        df2.columns = ['CHROM', 'START', 'END']        
        df2 = df2.astype({'START': 'int32'})
        df2 = df2.astype({'END': 'int32'})
        
        output_name = sample + '.cov_cut_mergedadjacent.bed'
        df2.to_csv(os.path.join(output, output_name), index = False, header = False, sep = '\t')


def return_files_in_dir(path, ext='/*.coverage'):
    return glob.glob(path + ext, recursive=True)

def import_mutect_vcf(path, skiprows_num):
    return pd.read_csv(path, sep='\t', skiprows=skiprows_num)

def add_filename_col(file, df):
    df['sampleId'] = ntpath.basename(file)
    return df

def import_vcf(file):
    with open(file) as f:
        content = f.readlines()
    df = pd.read_csv(file, sep='\t', skiprows=[i for i, s in enumerate(content) if '#CHROM\tPOS\tID\tREF\tALT\t' in s][0])
    return df

def import_append_process_vcfs(files):
    df = pd.DataFrame()
    for file in files:
        with open(file) as f:
            content = f.readlines()
        df_new = import_mutect_vcf(file, [i for i, s in enumerate(content) if '#CHROM\tPOS\tID\tREF\tALT\t' in s][0])
        df_new = add_filename_col(file, df_new)
        df = df.append(df_new)
    return df.reset_index(drop=True)

def save_vcf_header(file):
    with open(file) as f:
        content = f.readlines()
    info = [s for i, s in enumerate(content) if '#' in s]
    return info

def remove_indels(vcfs):
    return vcfs[(vcfs['REF'].apply(lambda x: len(x)) == 1) & (vcfs['ALT'].apply(lambda x: len(x)) == 1)].reset_index(drop=True)

def prep_annovar(vcfs, output_path):
    annot = vcfs.copy()
    csv_path = os.path.join(output_path, 'annot.csv')
    annot.to_csv(csv_path, header=False, sep='\t', index=False)
    
def import_merge_annovar_annotations(vcfs, path_intermed_in):
    annot = pd.read_csv(path_intermed_in)
    vcfs = vcfs.reset_index(drop=True)
    vcfs = pd.concat([vcfs, annot[annot.columns[5:]]], axis=1)
    return vcfs       

def gen_dummies(vcfs, Genotyped = True):
    # Annotate PASS Variants
    vcfs['PASS'] = vcfs['FILTER']== 'PASS'
    # Annotate Common Variants w/ MAF>1% (1000G, ExAC, ESP5600)
    vcfs['Common Variant'] = (vcfs['1000g2015aug_all'] > 0.01) | (vcfs['ExAC_ALL'] > 0.01) | (vcfs['esp6500siv2_all'] > 0.01)    # Annotate Single Cell Variants
    vcfs['1000G_blacklist'] = vcfs['bed4'].notnull()
    if Genotyped: 
        vcfs['CE_Indel'] = vcfs['bed6'].notnull()      
    return vcfs


if __name__ == "__main__":
    main()
