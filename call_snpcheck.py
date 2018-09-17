#!/usr/bin/env python

import os
import glob
import logging
import subprocess

from ngsscriptlibrary import parse_ngssnpcheck

HOME = os.path.expanduser('~')
JAVA = os.path.join(HOME, 'programs', 'jre1.8.0_144', 'bin', 'java')
GATK = os.path.join(HOME, 'programs', 'GenomeAnalysisTK-3.8-0-ge9d806836', 'GenomeAnalysisTK.jar')
REF = os.path.join(HOME, 'Documents', 'referentie', 'hg19.fa')
ALLELES = os.path.join(HOME, 'Documents', 'ngstargets', 'varia', 'NGS-SNPcheck.vcf')

def run_command(command):
    subprocess.call(command, shell=True, stderr=logging.INFO)

def call_snpcheck(bam):
    command = f'''{JAVA} -jar {GATK} -R {REF} -T HaplotypeCaller \
    -I {bam} -o {bam}.vcf -L {ALLELES} -alleles {ALLELES} \
    -gt_mode GENOTYPE_GIVEN_ALLELES -GQB 0
    '''
    try:
        run_command(command)
    except(OSError, CalledProcessError) as e:
        logging.error(f'{e}')
    else:
        snp_dict = parse_ngssnpcheck(f'{i}.vcf')

    return snp_dict

    
if __name__ == '__main__':

    logging.basicConfig(level=logging.DEBUG, filename='example.log')    
    
    run_command('map_amplicons.sh')
    
    for i in glob.glob(os.path.join(os.getcwd(), '*.bam')):
        logging.info(f'{i} verwerken')
        snpout = call_snpcheck(i)
        if not snpout:
            continue
            
        with open(f'{i}.ampliconngssnpcheck.txt', 'w') as f:
            for k, v in snpout.items():
              f.write(f'{k}\t{v}\n')       
    

