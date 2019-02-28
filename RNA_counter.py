import os
import argparse
import pybedtools


def get_exons_from_bed(bed_file):
    exons = list()
    with open(bed_file) as f:
        for line in f:
            if line in ['\n', '\r\n']:
                continue
            chrom, start, end, info, *_ = line.split()
            info = info.split('_')
            exon = info[3]
            exons.append(exon)
    return exons


def get_coverage_per_exon(bam_file, bed_file):
    bed = pybedtools.BedTool(bed_file)
    bam = pybedtools.BedTool(bam_file)

    coverage = bed.coverage(bam, split=True, counts=True)

    total = 0
    counts = list()
    
    
    for _ in coverage:
        chrom, start, end, info, count = _
        total += int(count)
        counts.append(int(count)) 
    
    counts_percentage = [count / total for count in counts]

    return counts_percentage


def parse_bam_files(bam_files, bed_file):
    output = dict()

    for bam_file in bam_files:
        bam_file_name = os.path.basename(bam_file)
        sample = bam_file_name.split('.')[0]
        output[sample] = get_coverage_per_exon(bam_file, bed_file)

    return output


def main(bam_files, bed_file, output_file):   
  
    exons = get_exons_from_bed(bed_file)    
    output = parse_bam_files(bam_files, bed_file)
    samples = sorted(output.keys())
    
    with open(output_file, 'w') as f:
        f.write('Exon\t')
        f.write('\t'.join(samples))
        f.write('\n')
        for i, exon in enumerate(exons):
            f.write(f'{exon}\t')
            [f.write(f'{output[sample][i]}\t') for sample in samples]
            f.write('\n')            


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--bamfiles", type=str, metavar='', nargs='*',
                        required=True, help="BAM files")
    parser.add_argument("--bedfile", type=str, metavar='',
                        required=True, help="DoC file(s) van controle(s)")
    parser.add_argument("-o", "--output", type=str, metavar='',
                        default='RNA_counts.txt',
                        help="Output file name")

    args = parser.parse_args()
    bam_files = args.bamfiles
    bed_file = args.bedfile
    output_file = args.output
    main(bam_files, bed_file, output_file)            

