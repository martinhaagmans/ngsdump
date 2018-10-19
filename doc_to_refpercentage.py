# coding: utf-8
import csv
import pysam
import argparse
import pandas as pd

from collections import namedtuple

def parse_doc(fn, ref, loci):
    data = dict()
    with open(fn) as f:
        spamreader = csv.reader(f, delimiter='\t')
        _header = next(spamreader)
        for line in spamreader:
            nonref = list()
            locus, _TD, _ADS, DP, basecounts = line
            if locus not in loci:
                continue
            refbase = ref[locus]
            for bases in basecounts.split(' '):
                base, cov = bases.split(':')
                cov = int(cov)
                try:
                    basep = int(cov) / int(DP)
                except ZeroDivisionError:
                    basep = 0
                if refbase == base:
                    refp = basep
                elif refbase != base:
                    nonref.append((base, basep))
            if int(DP) > 100:
                locus_data = namedtuple('mosaic_out', 'DP, refpercentage, nonreflist')
                data[locus] = locus_data(int(DP), refp, nonref)

    return data


# Add workdir to DoC file

def get_loci_from_docfile(docfile):
    loci = list()
    with open(docfile) as f:
        header = next(f)
        for line in f:
            locus, *_ = line.split()
            loci.append(locus)
    return loci


def get_ref_dict(REFERENCE):
    refd = dict()
    for locus in loci:
        base = pysam.faidx(REFERENCE,
                           '{}:{}-{}'.format(locus.split(':')[0],
                                             locus.split(':')[1],
                                             locus.split(':')[1])
                           ).split('\n')[1]
        refd[locus] = base
    return refd


def create_output(outfile, patient, control):
    with open(outfile, 'w') as f_out:
        for k, v in patient.items():
            try:
                c = control[k]
            except KeyError:
                c = 'GeenCoverage'
            f_out.write(f'{k}\t{v.DP}\t{100*v.refpercentage:.3f}\t{v.nonreflist}\t{c.DP}\t{100*c.refpercentage:3f}\t{c.nonreflist}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample", type=str, metavar='',
                        required=True, help="DoC file(s) van sample(s)")
    parser.add_argument("-c", "--control", type=str, metavar='',
                        required=True, help="DoC file(s) van controle(s)")
    parser.add_argument("-r", "--reference", type=str, metavar='',
                        default='/home/mahaagmans/Documents/referentie/hg19.fa',
                        help="Reference fasta")
    parser.add_argument("-o", "--output", type=str, metavar='',
                        default='final_output.txt',
                        help="Reference fasta")

    args = parser.parse_args()

    doc_sample = args.sample
    doc_control = args.control

    loci = get_loci_from_docfile(doc_sample)
    refd = get_ref_dict(args.reference)

    data_sample = parse_doc(args.sample, refd, loci)
    data_control = parse_doc(args.control, refd, loci)

    create_output(args.output, data_sample, data_control)
