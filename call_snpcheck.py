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
    subprocess.call(command, shell=True)


def call_snpcheck(bam):
    "Call SNP's with GATK, parse VCF and return dictionary"

    command = f"""{JAVA} -jar {GATK} -R {REF} -T HaplotypeCaller \
    -I {bam} -o {bam}.vcf -L {ALLELES} -alleles {ALLELES} \
    -gt_mode GENOTYPE_GIVEN_ALLELES -GQB 0
    """
    snp_dict = dict()
    
    try:
        run_command(command)
    except(OSError) as e:
        print(e)
    else:
        snp_dict = parse_ngssnpcheck(f'{i}.vcf')
    
    return snp_dict


if __name__ == '__main__':

   run_command('map_amplicons.sh')

    for i in glob.glob(os.path.join(os.getcwd(), '*.bam')):
        print(i)
        snpout = call_snpcheck(i)
        if not snpout:
            continue

        with open(f'{i}.ampliconngssnpcheck.txt', 'w') as f:
            for k, v in snpout.items():
                f.write(f'{k}\t{v}\n')
