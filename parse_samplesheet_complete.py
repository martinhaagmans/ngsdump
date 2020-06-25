import os
import csv 
import sys

ngslib = os.path.join("C:\\", "Users", "mahaa", "Documents", "GitHub", "ngsscriptlibrary")

sys.path.append(ngslib)

from ngsscriptlibrary import TargetDatabase


def get_header(samplesheet):
    "Read samplesheet and find line with Sample_ID. Return integer."
    with open(samplesheet, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('Sample_ID'):
                return i


def samplesheet_to_sample_genesis(samplesheet):
    """Read samplesheet part with sample-genesis info, return list of tuples."""
    samples = list()
    with open(samplesheet) as f:
        reader = csv.reader(f)
        for i, line in enumerate(reader):
            if i > get_header(samplesheet):
                samples.append(line)
    sample_genesis = list()
    for line in samples:
        sample = line[0]
        genesis = line[-1]
        genesis = genesis.replace('.NGS', '')
        genesis = genesis.replace('.SV', '')
        sample_genesis.append((sample, genesis))
    return sample_genesis


def parse_samplesheet_for_pipeline(samplesheet, db, exclude=None):
    """Read samplesheet part with sample info. Use ngsscriptlibrary to get
    targets to use and analyses to perform. Return dict.
    """
    samples = list()
    if exclude is None:
        exclude = list()

    with open(samplesheet) as f:
        reader = csv.reader(f)
        for i, line in enumerate(reader):
            if i > get_header(samplesheet):
                samples.append(line)
    TD = TargetDatabase(db)
    samples_todo = dict()
    for line in samples:
        genesis = line[-1]
        serie = line[-2]
        if genesis in exclude or serie in exclude:
            continue
        genesis = genesis.replace('.NGS', '')
        genesis = genesis.replace('.SV', '')
        samples_todo[line[0]] = TD.get_todo(genesis)
        samples_todo[line[0]]['serie'] = serie
        samples_todo[line[0]]['genesis'] = genesis
        if not samples_todo[line[0]]['amplicon']:
            vcapture = samples_todo[line[0]]['capture']
            oid = TD.get_oid_for_vcapture(vcapture)
            samples_todo[line[0]]['oid'] = oid
    return samples_todo



def get_file_locations(todo, targetrepo):
    """Read dict with targets and analyses and add correct file locations.
    Return dict
    """
    for s in todo.keys():
        picard = '{}_target.interval_list'.format(todo[s]['capture'])
        picard = os.path.join(targetrepo, 'captures', picard)

        cnvtarget = '{}_target.bed'.format(todo[s]['capture'])
        cnvtarget = os.path.join(targetrepo, 'captures', cnvtarget)

        cap_is_pakket = todo[s]['capispakket']

        if cap_is_pakket:
            pakket_name = todo[s]['capture']
            annot = '{}_target.annotated'.format(pakket_name)
            annot = os.path.join(targetrepo, 'captures', annot)
            varcal = '{}_generegions.bed'.format(pakket_name)
            varcal = os.path.join(targetrepo, 'captures', varcal)
            sanger = '{}_target.bed'.format(pakket_name)
            sanger = os.path.join(targetrepo, 'captures', sanger)
            pakket = '{}_target.bed'.format(pakket_name)
            pakket = os.path.join(targetrepo, 'captures', pakket)
        elif not cap_is_pakket:
            pakket_name = todo[s]['pakket']
            annot = '{}_target.annotated'.format(pakket_name)
            annot = os.path.join(targetrepo, 'pakketten', annot)
            varcal = '{}_generegions.bed'.format(pakket_name)
            varcal = os.path.join(targetrepo, 'pakketten', varcal)
            sanger = '{}_target.bed'.format(pakket_name)
            sanger = os.path.join(targetrepo, 'pakketten', sanger)
            pakket = '{}_target.bed'.format(pakket_name)
            pakket = os.path.join(targetrepo, 'pakketten', pakket)
        if todo[s]['panel'] is not None:
            annot = '{}_target.annotated'.format(todo[s]['panel'])
            annot = os.path.join(targetrepo, 'panels', annot)
            sanger = '{}_target.bed'.format(todo[s]['panel'])
            sanger = os.path.join(targetrepo, 'panels', sanger)
        if todo[s]['riskscore']:
            riskscorevcf = '{}_riskscore.vcf'.format(todo[s]['genesis'])
            riskscorevcf = os.path.join(targetrepo, 'varia', riskscorevcf)
            todo[s]['riskscorevcf'] = riskscorevcf

        todo[s]['annot'] = annot
        todo[s]['picard'] = picard
        todo[s]['cnvtarget'] = cnvtarget
        todo[s]['pakkettarget'] = pakket
        todo[s]['varcal'] = varcal
        todo[s]['sanger'] = sanger

    return todo
