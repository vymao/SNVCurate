import argparse
import os
import re
import ntpath
import sys
import glob

def parse_args():
    """Uses argparse to enable user to customize script functionality""" 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-mutect_path', default=None, help='path to mutect output')
    parser.add_argument('-muse_path', default=None, help='path to muse output')
    parser.add_argument('-haplotypecaller_path', default=None, help='path to haplotypecaller output')
    parser.add_argument('-bam_path', default=None, help='path to BAM output')
    parser.add_argument('-csv', default=None, help='directory to which the output directory "/.Mutect2/" will be written to')
    parser.add_argument('-out', default='./')

    return parser.parse_args()


def main(): 
    dirname = os.path.dirname(os.path.abspath(__file__))

    args = vars(parse_args())
    if args['out'] == '.' or args['out'] == './': args['out'] = os.getcwd() 

    if args['csv'] is None: 
        print('No csv given.')
        sys.exit()
    if args['mutect_path'] is None: 
        print('No Mutect2 output path given.')
        sys.exit()

    tumor_normal_matched = os.path.join(args['out'], 'Tumor_Normal.csv')

    with open(args['csv'], 'r') as main: 
        for index, line in enumerate(main): 
            if index == 0: 
                continue
            else: 
                file_list = line.rstrip().split(',')
                sample = line[0].rstrip().split('.')[0]

                os.symlink(line[0], os.path.join(args['bam_path'], sample + '.bam'))
                matched = sample + '.bam'

                os.makedirs(os.path.join(args['mutect_path'], sample), exist_ok=True)
                mutect_sample_path = os.path.join(args['mutect_path'], sample)
                os.symlink(line[1], os.path.join(mutect_sample_path, sample + '.vcf'))

                if args['muse_path'] is not None: 
                    if line[2] is not None: 
                        os.makedirs(os.path.join(args['muse_path'], sample), exist_ok=True)
                        muse_sample_path = os.path.join(args['muse_path'], sample)
                        os.symlink(line[2], os.path.join(muse_sample_path, sample + '.vcf'))
                if args['haplotypecaller_path'] is not None: 
                    if line[3] is not None: 
                        normal_sample = line[3].rstrip().split('.')[0]
                        os.symlink(line[3], os.path.join(args['haplotypecaller_path'], normal_sample + '.vcf'))

                        matched += ',' + normal_sample + '.bam'
                    else: 
                        matched += ',N/A'
                else: 
                    matched += ',N/A'
                with open(tumor_normal_matched, 'a') as csv_out: 
                    csv_out.write(matched)


if __name__ == "__main__":
    main()

