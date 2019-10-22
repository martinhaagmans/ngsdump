# coding: utf-8
import os
import json
import logging
import argparse
import pybedtools
from ngsscriptlibrary import TargetAnnotation
from ngsscriptlibrary import get_picard_header


def check_bed_intervals(bed):
    bed_intervals_correct = True
    with open(bed) as f:
        for line in f:
            _c, s, e, *_ = line.split()
            if not s < e:
                bed_intervals_correct = False
    return bed_intervals_correct


def correct_intervals(bed):
    filepath, ext = os.path.splitext(bed)
    corrected_bed = '{}.corrected_intervals{}'.format(filepath, ext)
    with open(bed) as f, open(corrected_bed, 'w') as fout:
        for line in f:
            c, s, e, *_ = line.split()
            if s < e:
                fout.write('{}\t{}\t{}\n'.format(c, s, e))
            elif s > e:
                fout.write('{}\t{}\t{}\n'.format(c, e, s))
            elif s == e:
                fout.write('{}\t{}\t{}\n'.format(c, int(s) - 1, e))
    return corrected_bed


def merge_bed(bed_location):
    filepath, ext = os.path.splitext(bed_location)
    merged_bed_location = ('{}.merged{}'.format(filepath, ext))

    if not check_bed_intervals(bed_location):
        bed_location = correct_intervals(bed_location)

    bed = pybedtools.BedTool(bed_location)
    bed = bed.sort()
    bed = bed.merge()
    bed.saveas(merged_bed_location)
    return merged_bed_location


def check_if_all_intervals_annotated(bed):
    not_annotated_intervals = list()
    for line in bed:
        try:
            _c, _s, _e, _annotation, *_ = line
        except ValueError:
            not_annotated_intervals.append(line)
    return not_annotated_intervals


def get_genes_in_bed(bed):
    genes_in_bed = list()
    for line in bed:
        _c, _s, _e, gene, *_ = line
        genes_in_bed.append(gene)
    return genes_in_bed


def get_genes_from_request_not_in_bed(genes_in_bed, genes_requested):
    not_in_bed = list()
    for gene in genes_requested:
        if gene not in genes_in_bed and gene not in not_in_bed:
            not_in_bed.append(gene)
    return not_in_bed


def get_genes_from_bed_not_in_request(genes_in_bed, genes_requested):
    not_in_request = list()
    for gene in genes_in_bed:
        if gene not in genes_requested and gene not in not_in_request:
            not_in_request.append(gene)
    return not_in_request


def get_nogene_targets(bed):
    nogene_regions = list()
    for line in bed:
        _c, _s, _e, gene, *_ = line
        if gene == 'NOGENE':
            nogene_regions.append(line)
    return nogene_regions


def write_report(genes_not_in_bed, genes_not_in_list, nogene_targets, targetname):

    with open('{}_report.txt'.format(targetname), 'w') as f:
        f.write('Wel in lijst, niet in BED:\n')
        if len(genes_not_in_bed) > 0:
            f.write('{}\n\n'.format(', '.join(genes_not_in_bed)))
        else:
            f.write('Geen\n\n')

        f.write('Wel in BED, niet in lijst:\n')
        if len(genes_not_in_list) > 0:
            f.write('{}\n\n'.format(', '.join(genes_not_in_list)))
        else:
            f.write('Geen\n\n')

        f.write('Niet in genregio:\n')
        if len(nogene_targets) > 0:
            for target in nogene_targets:
                chromosome, start, end, _gene = target
                f.write('{}:{}-{}\n'.format(chromosome, start, end))
        else:
            f.write('Geen\n\n')


def write_annotated_bed(bed, targetname):
    with open('{}.annotated'.format(targetname), 'w') as f:
        for line in bed:
            chromosome, start, end, gene = line
            f.write('{}\t{}\t{}\t{}\n'.format(chromosome, start, end, gene))
    

def write_generegion_bed(genes, targetname, nogene_regions):
    with open('{}_generegions.bed'.format(targetname), 'w') as f:
        regions_done = list()
        for gene in genes:
            for region in TA.get_region(gene):
                if region in regions_done:
                    continue
                regions_done.append(region)
                chromosome, start, end = region
                f.write('{}\t{}\t{}\t{}\n'.format(chromosome, start, end, gene))
        for region in nogene_regions:
            chromosome, start, end, gene = region
            f.write('{}\t{}\t{}\t{}\n'.format(chromosome, start, end, gene))


def write_picard_targetfile(bed, targetname):
    picard_header = get_picard_header()
    with open('{}_target.interval_list'.format(targetname), 'w') as f:
        for line in picard_header:
            f.write('{}'.format(line))
        for line in bed_annotated:
            chromosome, start, end, _gene = line
            f.write('{}\t{}\t{}\t+\t{}\n'.format(chromosome, start, end, targetname))


def write_genelist_to_json(genes, targetname):
    with open('{}_genes.json'.format(targetname), 'w') as f:
        json.dump(genes, f)


if __name__ == "__main__":
    basedir = os.path.dirname(os.path.realpath(__file__))
    log = os.path.join(basedir, 'samplesheets.log')
    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    parser = argparse.ArgumentParser(description="Check BED for targetrepo")
    parser.add_argument('-b', '--bedfile', type=str, required=True,
                        help='BED-file to check')
    parser.add_argument('-g', '--genes', type=str,
                        help='Genelist to check')

    args = parser.parse_args()
    genelist_inputfile = args.genes
    bed_input = args.bedfile
    bed_merged = merge_bed(bed_input)
    logger.info('BED merged.')

    base = os.path.basename(bed_input)
    targetname, ext = os.path.splitext(base)

    if genelist_inputfile:
        TA = TargetAnnotation(bed_merged, genes=genelist_inputfile)
        bed_annotated = TA.annotate_bed_and_filter_genes()

        genelist = TA.genes
        nogene_targets = get_nogene_targets(bed_annotated)
        genes_in_bed = get_genes_in_bed(bed_annotated)
        genes_not_in_bed = get_genes_from_request_not_in_bed(genes_in_bed, genelist)
        genes_not_in_list = get_genes_from_bed_not_in_request(genes_in_bed, genelist)

        write_report(genes_not_in_bed, genes_not_in_list, nogene_targets, targetname)
        write_generegion_bed(genelist, targetname, nogene_targets)
        logger.info('Generegions BED created.')
        write_genelist_to_json(genelist, targetname)

    else:
        TA = TargetAnnotation(bed_merged, genes=genelist_inputfile)
        bed_annotated = TA.annotate_bed_and_filter_genes()

    write_annotated_bed(bed_annotated, targetname)
    logger.info('Annotated BED created.')
    write_picard_targetfile(bed_annotated, targetname)
    logger.info('Picard target created.')
