# coding: utf-8
import os
import csv
import pysam
import argparse
import pandas as pd

from collections import namedtuple

def remove_extension(fn):
    fnbase, fnext = os.path.splitext(fn)
    return fnbase

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

            locus_data = namedtuple('mosaic_out', 'DP, refpercentage, nonreflist')
            data[locus] = locus_data(int(DP), refp, nonref)

    return data


# Add workdir to DoC file

def get_loci_from_docfile(docfile):
    loci = list()
    with open(docfile) as f:
        _header = next(f)
        for line in f:
            locus, *_ = line.split()
            loci.append(locus)
    return loci


def get_ref_dict(ref, loci):
    refd = dict()
    for locus in loci:
        base = pysam.faidx(ref,
                           '{}:{}-{}'.format(locus.split(':')[0],
                                             locus.split(':')[1],
                                             locus.split(':')[1])
                           ).split('\n')[1]
        refd[locus] = base
    return refd


def create_output(outfile, sample_data, control_data, loci, allIDs):
    allout = dict()

    for locus in loci:

        out = list()

        for data in sample_data:
            out.append(data[locus])

        for data in control_data:
            out.append(data[locus])

        allout[locus] = out

    with open(outfile, 'w') as f_out:
        f_out.write('Locus')

        for id in allIDs:
            f_out.write(f'\t{id} DP\t{id} refP\t{id} nonref')
        f_out.write('\n')

        for locus in loci:
            f_out.write(f'{locus}')

            for _ in allout[locus]:
                f_out.write(f'\t{_.DP}\t{100*_.refpercentage:.3f}\t{_.nonreflist}')

            f_out.write('\n')

def collect_data(doc_sample, doc_control, loci, ref):
    
    refd = get_ref_dict(ref, loci)

    data_sample = list()
    data_control = list()
    allIDs = list()

    if len(doc_sample) == 1:
        data_sample.append(parse_doc(doc_sample[0], refd, loci))
        allIDs.append(remove_extension(doc_sample))
    else:
        for doc in doc_sample:
            data_sample.append(parse_doc(doc, refd, loci))
            allIDs.append(remove_extension(doc))

    if len(doc_control) == 1:
        data_control.append(parse_doc(doc_control[0], refd, loci))
        allIDs.append(remove_extension(doc_control[0]))
    else:
        for doc in doc_control:
            data_control.append(parse_doc(doc, refd, loci))
            allIDs.append(remove_extension(doc))
    return data_sample, data_control, allIDs


def main(doc_sample, doc_control, ref, output):
    loci = get_loci_from_docfile(doc_sample[0])
    data_sample, data_control, allIDs = collect_data(doc_sample, doc_control, loci, ref)
    create_output(output, data_sample, data_control, loci, allIDs)    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample", type=str, metavar='', nargs='*',
                        required=True, help="DoC file(s) van sample(s)")
    parser.add_argument("-c", "--control", type=str, metavar='', nargs='*',
                        required=True, help="DoC file(s) van controle(s)")
    parser.add_argument("-r", "--reference", type=str, metavar='',
                        default='/home/mahaagmans/Documents/referentie/hg19.fa',
                        help="Reference fasta")
    parser.add_argument("-o", "--output", type=str, metavar='',
                        default='final_output.txt',
                        help="Output file name")

    args = parser.parse_args()
    output = args.output
    ref = args.reference
    doc_sample = args.sample
    doc_control = args.control

    main(doc_sample, doc_control, ref, output)

