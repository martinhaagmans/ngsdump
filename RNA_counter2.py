import glob
import argparse
import subprocess

import pysam
import pybedtools

from math import floor
from statistics import mean, stdev


def get_reads(sample, readdir):
    r1 = glob.glob(f'{readdir}/{sample}_*R1*.gz')[0]
    r2 = glob.glob(f'{readdir}/{sample}_*R2*.gz')[0]
    return r1, r2


def align_reads(sample, r1, r2, ref, amplicon=True):
    bamfile = f'{sample}.sorted.bam'
    command = f"""(hisat2 -x {ref} -1 {r1} -2 {r2} --threads 2 |\
    samtools view -hu - |\
    samtools sort -o {bamfile}) > {sample}.hisat2.log 2>&1 \
    && samtools index {bamfile}
    """
    subprocess.run(command, shell=True)
    return bamfile
    

def usable_read(read):
    "Check read and return boolean"
    if read.is_duplicate:
        return False
    elif read.is_unmapped:
        return False
    elif read.is_supplementary:
        return False
    else:
        return True


def get_reads_per_locus(bamfile):
    reads_per_locus = dict()

    alignments = pysam.AlignmentFile(bamfile, "rb")

    for read in alignments:
        
        if not usable_read(read):
            continue

        chromosome = read.reference_name
        read_name = read.query_name
        for position in read.get_reference_positions():
            locus = f'{chromosome}:{position}'
            if locus in reads_per_locus:
                if not read_name in reads_per_locus[locus]:
                    reads_per_locus[locus].add(read_name)
            else:
                reads_per_locus[locus] = {read_name}
    alignments.close()

    print(f'{bamfile} done.')
    return reads_per_locus


def parse_reads_per_locus(reads_per_locus):
    bedlist = list()
    for locus, reads in reads_per_locus.items():
        count = len(reads)
        if int(count) < 10:
            continue
        chromosome, position = locus.split(':')
        position = int(position)
        bedlist.append((chromosome, position - 1, position, count))

    bed = pybedtools.BedTool(sorted(bedlist))
    bed.sort()
    bed_merged = bed.merge()
    bed_intersect = bed_merged.intersect(bed, wa=True, wb=True)

    interval_average = list()

    for interval in bed_merged:
        counts = list()
        interval = str(interval).rsplit()
        interval_chromosome, interval_start, interval_end = interval
        for line in bed_intersect:
            line = str(line).rsplit()
            chromosome, start, end, *_ = line
            
            # get booleans
            chrom = (interval_chromosome == chromosome)
            start = (interval_start == start)
            end = (interval_end == end)

            if chrom and start and end:
                count = line[-1]
                counts.append(int(count))

        count_mean = floor(mean(counts))
        count_stddev = floor(stdev(counts))
        interval_average.append((interval_chromosome, interval_start, interval_end, count_mean, count_stddev))

    bed_counts = pybedtools.BedTool(sorted(interval_average))
    return bed_counts

    
def get_counts_per_region(ucsc_bed, bamfile):
    reads_per_locus = get_reads_per_locus(bamfile)
    bed = parse_reads_per_locus(reads_per_locus)
    bed_uscs = pybedtools.BedTool(ucsc_bed)
    return bed_uscs.intersect(bed, wao=True)


def write_bed_counts_per_region(sample, bed_counts_per_region):
    with open(f'{sample}_count.bed', 'w') as f:
        for line in bed_counts_per_region:
            line = str(line).rsplit()

            (ucsc_chrom, ucsc_start, ucsc_end, nmnr, exon, 
              _, _, _ ,count, count_stdev, overlap) = line
            if count == '.':
                count = 0
            else:
                count = int(count)
            if count_stdev == '-1':
                count_stdev = 0
            else:
                count_stdev = int(count_stdev)

            f.write(f'{ucsc_chrom}\t{ucsc_start}\t{ucsc_end}\t'
                    f'{nmnr}\t{exon}\t{count}\t{count_stdev}\n')    


if __name__ == '__main__':

    REF = "/home/mahaagmans/Documents/referentie/hisat/grch37_tran/genome_tran"
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True, help="SampleID(s)")
    parser.add_argument("--reads", required=True, help="Directory with read files for sample")
    parser.add_argument("--bedfile", required=True, help="UCSC bed")
    parser.add_argument("--reference", default=REF, help="Hisat reference")                        

    args = parser.parse_args()
    sample = args.sample
    reference = args.reference
    readdir = args.reads
    ucsc_bed = args.bedfile

    r1, r2 = get_reads(sample, readdir)
    bamfile = align_reads(sample, r1, r2, reference)
    
    counts_per_region = get_counts_per_region(ucsc_bed, bamfile)
    write_bed_counts_per_region(sample, counts_per_region)

