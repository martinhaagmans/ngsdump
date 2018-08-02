# coding: utf-8
import pysam
import pandas as pd

from ngsscriptlibrary.mosaic import Mosaic
from ngsscriptlibrary.mosaic import parse_doc
from ngsscriptlibrary.mosaic import bedfile_to_locilist

REFERENCE = '/home/mahaagmans/Documents/referentie/hg19.fa'
PATIENTDOC = '18D6648.DoC'
CONTROLDOC = 'NC1.DoC'

def get_loci_from_docfile(docfile):
    loci = list()
    with open(docfile) as f:
        header = next(f)
        for line in f:
            locus, *_ = line.split('\t')
            loci.append(locus)
        

loci = get_loci_from_docfile(PATIENTDOC)

refd = dict()

for locus in loci:
    base = pysam.faidx(REFERENCE,
                       '{}:{}-{}'.format(locus.split(':')[0],
                                         locus.split(':')[1],
                                         locus.split(':')[1])
                       ).split('\n')[1]
    refd[i] = base
                 
patient = dict()
control = dict()

for k, v in parse_doc(PATIENTDOC, refd, loci).items():
    if v.DP > 100 :
        patient[k] = v
        
      
for k, v in parse_doc(CONTROLDOC, refd, loci).items():
    if v.DP > 100 :
        control[k] = v
        

with open('final_output.txt', 'w') as f_out:
    for k, v in patient.items():
        try:
            c = control[k]
        except KeyError:
            c = 'GeenCoverage'
        f_out.write(f'{k}\t{v.DP}\t{100*v.refpercentage:.3f}\t{v.nonreflist}\t{c.DP}\t{100*c.refpercentage:3f}\t{c.nonreflist}\n')

