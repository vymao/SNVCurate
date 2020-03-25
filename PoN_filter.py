import os
import pandas as pd
import ntpath
import glob
import re
import warnings
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline
import subprocess
import argparse
import sys
import pysam

def parse_args():
    """Uses argparse to enable user to customize script functionality"""
    parser = argparse.ArgumentParser()
    parser.add_argument('-somatic_vcf', help='path to initial somatic SNVs from initial filtering')
    parser.add_argument('-normal_vcf', default=None, help='path to normal file for InDel mask creation')
    parser.add_argument('-annovar', help='path to annovar database')
    parser.add_argument('-reference', default='hg19')
    parser.add_argument('-pon', default='/n/data1/hms/dbmi/park/victor/references/TCGA_1000_PON.hg19.REORDERED.vcf')
    parser.add_argument('-bam')
    parser.add_argument('-out')


    return parser.parse_args()

def main():
    args = vars(parse_args())
    path_vcfs_intersection = os.path.dirname(args['somatic_vcf'])
    basename = os.path.basename(args['somatic_vcf']).split('.')[0]


    path_intermed_out = args['annovar']
    somatic_file = args['somatic_vcf']

    vcfs = import_append_process_vcfs(somatic_file)
    vcfs = vcfs[vcfs['#CHROM'].apply(lambda x: len(str(x)) <=5)]
    vcfs = vcfs[vcfs['FILTER']=='PASS']
	prep_annovar(vcfs, basename + '.annot_normal.csv', path_vcfs_intersection)

	if args['normal_vcf'] is not None: 
		path_vcfs_germline = os.path.dirname(args['normal_vcf'])
		germline = import_append_process_vcfs(args['normal_vcf'])
		germline_indel_vcf = keep_only_indels(germline)
		create_indel_mask(germline_indel_vcf, os.path.join(args["annovar"], 'FCD_indels.bed'))
		output_name = os.path.join(path_vcfs_intersection, basename + "_Common_filtered")
		sample_paths = os.path.join(path_vcfs_intersection, basename+ '.annot_normal.csv')

		genotyped = True

		command = "perl /home/mk446/bin/annovar/table_annovar.pl " + sample_paths + " " + args["annovar"] + " -buildver " + args["reference"] + " -out " + \
					output_name + " -remove -protocol refGene,1000g2015aug_all,exac03,esp6500siv2_all,bed,bed,bed,bed,bed,bed -operation g,f,f,f,r,r,r,r,r,r" + 
					" -bedfile simpleRepeat.bed,hg19_rmsk.bed,all_repeats.b37.bed,20141020.strict_mask.whole_genome.bed,all.repeatmasker.b37.bed,FCD_indels.bed " + 
					"--argument ',,,,-colsWanted 4,-colsWanted 5,,,,' -csvout -polish"
	else: 
		sample_paths = os.path.join(path_vcfs_intersection, basename+ '.annot_normal.csv')
		output_name = os.path.join(path_vcfs_intersection, basename + "_Common_filtered")

		genotyped = False

		command = "perl /home/mk446/bin/annovar/table_annovar.pl " + sample_paths + " " + args["annovar"] + " -buildver " + args["reference"] + " -out " + \
					output_name + " -remove -protocol refGene,1000g2015aug_all,exac03,esp6500siv2_all,bed,bed,bed,bed,bed -operation g,f,f,f,r,r,r,r,r" + 
					" -bedfile simpleRepeat.bed,hg19_rmsk.bed,all_repeats.b37.bed,20141020.strict_mask.whole_genome.bed,all.repeatmasker.b37.bed " + 
					"--argument ',,,,-colsWanted 4,-colsWanted 5,,,,' -csvout -polish"

	os.system(command)

	path_intermed_in = os.path.join(path_vcfs_intersection, output_name + '.hg19_multianno.csv')
	vcfs = import_merge_annovar_annotations(vcfs, path_intermed_in)
	vcfs = gen_dummies(vcfs, Genotyped)

	pon = import_vcf(path_pon)
	pon['PON'] = True
	pon = pon.drop_duplicates(['#CHROM', 'POS'])
	vcfs_merged = vcfs.merge(pon[['#CHROM', 'POS', 'PON']], how='left', on=['#CHROM', 'POS'])
	vcfs = vcfs_merged

	if genotyped:
    	vcfs = vcfs[(vcfs['Common Variant']==False) & (vcfs['1000G_blacklist']==True) & (vcfs['CE_Indel']==False) & (vcfs['PON']!=True)]
	else:
    	vcfs = vcfs[(vcfs['Common Variant']==False) & (vcfs['1000G_blacklist']==True) & (vcfs['PON']!=True)]

    path_bedtools_intersect = '/n/data1/hms/dbmi/park/alon/command_line_tools/Bedtools/Intersect.py'
    path_exome_capture_bed = '/n/data1/hms/dbmi/park/alon/files/FCD/Reference_Files/TruSeq_Exome_Targeted_Regions_Manifest/truseq-exome-targeted-regions-manifest-v1-2.bed'
    path_bins = '/n/data1/hms/dbmi/park/alon/files/FCD/Reference_Files/bins/bins.bed'
    path_bins_out = '/n/data1/hms/dbmi/park/alon/files/FCD/Filtration/SV/CNV_SegDup_Bins/bins_no_contigs.bed'
	path_bams = args['bam']
	path_soft_clipped_cutoff_out = os.path.join(path_vcfs_intersection, 'soft_clipped_cutoff.csv')

	if args['reference'] == 'hg19': hg19 = True
	else: hg19 = False


	clean_bins(path_bins, path_bins_out, genotyped)
	bins = pd.read_csv(path_bins_out, sep='\t', header=None)

	soft_clipped_cutoff = generate_soft_clipped_cutoff(path_vcfs_intersection, bins)
	soft_clipped_cutoff.to_csv(path_soft_clipped_cutoff_out, index=False)

	vcfs = merge_soft_clipped_cutoff(vcfs, soft_clipped_cutoff)
	vcfs = annotate_clipped_reads(vcfs, hg19)
	vcfs = vcfs[vcfs['clipped_reads'] < vcfs['Soft Clipped Cutoff']]
	vcfs.to_csv(os.path.join(path_vcfs_intersection, 'Final_Callset.txt'), sep = '\t', index=False)

	output_name = output_name = os.path.join(path_vcfs_intersection, basename + ".Final_Callset")
	prep_annovar(vcfs, basename + '.annot_normal.csv', path_vcfs_intersection)
	command = "perl /home/mk446/bin/annovar/table_annovar.pl " + os.path.join(path_vcfs_intersection, basename + '.annot_normal.csv') + " " + "/home/mk446/bin/annovar/humandb/" 
			+ " -buildver " + args["reference"] + " -out " + output_name + " -remove -protocol refGene,clinvar_20190305,dbnsfp33a -operation g,f,f -nastring . -vcfinput -polish"

	os.system(command)





def return_files_in_dir(path, ext='*somatic_variants_filtered_1.vcf'):
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

def import_append_process_vcfs(file):
    df = pd.DataFrame()
    with open(file) as f:
        content = f.readlines()
    df_new = import_mutect_vcf(file, [i for i, s in enumerate(content) if '#CHROM\tPOS\tID\tREF\tALT\t' in s][0])
    df_new = add_filename_col(file, df_new)

    return df_new

def save_vcf_header(file):
    with open(file) as f:
        content = f.readlines()
    info = [s for i, s in enumerate(content) if '#' in s]
    return info

def remove_indels(vcfs):
    return vcfs[(vcfs['REF'].apply(lambda x: len(x)) == 1) & (vcfs['ALT'].apply(lambda x: len(x)) == 1)].reset_index(drop=True)

def prep_annovar(vcfs, name, path_vcfs_intersection):
    annot = vcfs.copy()
    annot['END'] = annot['POS']
    annot = annot[list(annot.columns[:2]) + ['END'] + list(annot.columns[3:5])]
    csv_path = os.path.join(path_vcfs_intersection, name)
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

def keep_only_indels(vcfs):
    return vcfs[(vcfs['REF'].apply(lambda x: len(x)) > 1) | (vcfs['ALT'].apply(lambda x: len(x)) > 1)].reset_index(drop=True)

def create_indel_mask(indels, path_bed_out):
    indels['END_POS'] = indels['POS'].astype(int) + pd.DataFrame([indels['REF'].apply(lambda x: len(x)), indels['ALT'].apply(lambda x: len(x))]).max()
    indels[['POS', 'END_POS']] = indels[['POS', 'END_POS']].astype(int)
    indels['START_POS'] = indels['POS'] - 200
    indels['END_POS'] = indels['END_POS'] + 200
    indels['indic'] = 'INDEL'
    indels[['#CHROM', 'START_POS', 'END_POS', 'indic']].to_csv(path_bed_out, sep='\t', index=False, header=None)
    #return test

def clean_bins(path_bins, path_bins_out, Genotyped):
    # Import bins
    bins = pd.read_csv(path_bins, sep='\t', header=None)
    # Convert to 1-based coordinate system (for SAM) from 0-based (BED)
    bins[1] = bins[1] + 1
    # Remove non canonical Chroms 1-22 + X,Y,MT
    bins = bins[bins[0].isin([i for i in bins[0].unique() if len(i) <= 5])]
    # Convert chromosome naming convention to comply with samtools
    if Genotyped:
        bins[0] = bins[0].apply(lambda x: x.strip('chr'))
        bins[0] = bins[0].apply(lambda x: re.sub('M', 'MT', x))
    bins.to_csv(path_bins_out, header=None, index=False, sep='\t')
    
def generate_soft_clipped_cutoff(path_vcfs_intersection, bins):
    bam_files = []
    soft_clipped = []
    """
    for file in return_files_in_dir(path_vcfs_intersection):
        path_bam = path_bams + ntpath.basename(file).split('_v_')[0]
        bam_files.append(ntpath.basename(path_bam))
        soft_clipped.append(return_soft_clipped_percentile(path_bam, bins))
    """
    
    bam_files.append(ntpath.basename(path_bams))
    soft_clipped.append(return_soft_clipped_percentile(path_bams, bins))
    soft_clipped_cutoff = pd.DataFrame({'sampleId': bam_files, 'Soft Clipped Cutoff': soft_clipped})
    return soft_clipped_cutoff

def return_soft_clipped_percentile(path, bins):
    samfile = pysam.AlignmentFile(path, "rb")
    sampled_bins = bins.sample(n=50000)
    soft_clipped = []
    lst = zip(list(sampled_bins[0].astype(str)), list(sampled_bins[1]), list(sampled_bins[2]))
    for chrom, start, stop in lst:
        num_reads = 0
        num_sc_reads = 0
        for read in samfile.fetch(chrom, start, stop):
            num_reads += 1
            lst = re.findall(r"[^\W\d_]+|\d+", str(read).split('\t')[5])
            if 'S' in lst:
                num_sc_reads += 1
        if num_reads > 0:
            soft_clipped.append(num_sc_reads/(num_reads))
    return np.percentile(np.array(soft_clipped), 95)

def return_num_clipped_reads(chrom, start, stop, path, sam_file, hg19):
    #def strip_chrom(chromosome):
    #    chrom_number = chromosome.split('chr')[1]
    #    return chrom_number
    #chrom = chrom.apply(strip_chrom)
    if not hg19:
        chrom = chrom.split('chr')[1]
    samfile = sam_file
    num_reads = 0
    num_sc_reads = 0
    
    for read in samfile.fetch(chrom, start, stop):
        num_reads += 1
        lst = re.findall(r"[^\W\d_]+|\d+", str(read).split('\t')[5])
        if 'S' in lst:
            num_sc_reads += 1
    if num_reads > 0:
        return num_sc_reads/num_reads
    else:
        return 0
    return num_sc_reads

def annotate_clipped_reads(vcfs, hg19):
    df = vcfs.copy()
    df['POS'] = df['POS'].astype(int)
    df['POS_END'] = df['POS'] + 100
    df['POS'] = df['POS'] - 99
    df[['POS', 'POS_END']] = df[['POS', 'POS_END']].astype(int)
    last_sampleId = ''
    sam_file = ''
    for index, row in df.iterrows():
        if row['sampleId'] != last_sampleId:
            last_sampleId = row['sampleId']
            #sam_file = pysam.AlignmentFile(path_bams + row['sampleId'].split('_v_')[0], "rb")
            sam_file = pysam.AlignmentFile(path_bams, "rb")
            df.at[index, 'clipped_reads'] = return_num_clipped_reads(str(row['#CHROM']), int(row['POS']), int(row['POS_END']), path_bams + row['sampleId'], sam_file, hg19)
        else:
            df.at[index, 'clipped_reads'] = return_num_clipped_reads(str(row['#CHROM']), int(row['POS']), int(row['POS_END']), path_bams + row['sampleId'], sam_file, hg19)        
    vcfs['clipped_reads'] = df['clipped_reads']
    return vcfs

def merge_soft_clipped_cutoff(vcfs, soft_clipped_cutoff, path_vcfs_intersection):
    vcfs['sampleId'] = vcfs['sampleId'].astype(str)
    soft_clipped_cutoff['sampleId'] = [ntpath.basename(file) for file in return_files_in_dir(path_vcfs_intersection)]
    soft_clipped_cutoff[['sampleId']] = soft_clipped_cutoff[['sampleId']].astype(str)
    vcfs = vcfs.merge(soft_clipped_cutoff, how='left', on='sampleId')
    return vcfs

if __name__ == "__main__":
    main()