import os
import argparse
import pybedtools


def correct_chromosome_names(bed_file):
    file_location, file_name = os.path.split(bed_file)
    file_name, file_extension = os.path.splitext(file_name)
    corrected_file = os.path.join(file_location, 
                                  f'{file_name}_corrected{file_extension}')
    correct = False
    with open(bed_file) as f_in, open(corrected_file, 'w') as f_out:
        for line in f_in:
            if not line:
                continue
            elif line.startswith('chr'):
                correct = True
                f_out.write(line.replace('chr', '', 1))
    if not correct:
        os.remove(corrected_file)
        corrected_file = bed_file
        
    return corrected_file
    


def get_regions_from_bed(bed_file):
    regions = list()
    with open(bed_file) as f:
        for line in f:
            if line in ['\n', '\r\n']:
                continue
            chrom, start, end, info, *_ = line.split()
            info = info.split('_')
            nm = f'NM_{info[1]}'
            region_type = info[2]
            region_number = info[3]
            regions.append((chrom, start, end, nm, region_number))
    return region_type, regions


def get_coverage_per_exon(bam_file, bed_file):
    bed = pybedtools.BedTool(bed_file)
    bam = pybedtools.BedTool(bam_file)

    coverage = bed.coverage(bam, split=True, counts=True)

    total = 0
    counts = list()
    current_nm = str()
    counts_percentage = list()

    for line in coverage:
        chrom, start, end, info, *_ = line
        info = info.split('_')
        count = line[-1]
        nm = f'NM_{info[1]}'

        if not current_nm:
            current_nm = nm
            total += int(count)
            counts.append(int(count))
        elif nm == current_nm:
            total += int(count)
            counts.append(int(count)) 
            
        elif nm != current_nm:
            counts_percentage = counts_percentage + [c / total for c in counts]
            counts = list()
            total = 0
            counts.append(int(count)) 
            total += int(count)
            current_nm = nm
    counts_percentage = counts_percentage + [c / total for c in counts]
    return counts_percentage


def parse_bam_files(bam_files, bed_file):
    output = dict()

    for bam_file in bam_files:
        bam_file_name = os.path.basename(bam_file)
        sample = bam_file_name.split('.')[0]
        output[sample] = get_coverage_per_exon(bam_file, bed_file)

    return output


def main(bam_files, bed_file, output_file, correct):   
    if correct:
        bed_file = correct_chromosome_names(bed_file)
    region_type, regions = get_regions_from_bed(bed_file)
    region_type = region_type.capitalize()  
    output = parse_bam_files(bam_files, bed_file)
    samples = sorted(output.keys())

    
    with open(output_file, 'w') as f:
        f.write(f'Chr\tStart\tEnd\tNM\t{region_type}\t')
        f.write('\t'.join(samples))
        f.write('\n')
        for i, region in enumerate(regions):
            chrom, start, end, nm, exon = region
            f.write(f'{chrom}\t{start}\t{end}\t{nm}\t{exon}\t')
            [f.write(f'{output[sample][i]}\t') for sample in samples]
            f.write('\n')            


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--bamfiles", type=str, nargs='*',
                        required=True, help="BAM files")
    parser.add_argument("--bedfile", type=str,
                        required=True, help="DoC file(s) van controle(s)")
    parser.add_argument("-o", "--output", type=str, metavar='',
                        default='RNA_counts.txt',
                        help="Output file name")
    parser.add_argument("--correct", action='store_true', default=True,
                        help="Correct chromosome names")

    args = parser.parse_args()
    bam_files = args.bamfiles
    bed_file = args.bedfile
    output_file = args.output
    correct = args.correct
    main(bam_files, bed_file, output_file, correct)            

