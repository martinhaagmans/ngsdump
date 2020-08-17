import os
import sys
import json
import logging
import sqlite3
import argparse

logging.basicConfig(format='%(levelname)s: %(message)s')



def find_snpcheck(sample_to_check, serie_to_check, db):
    alt_loci = ['chr15:63351840', 'chr10:18789724', 'chr7:91726927', 'chr3:53700550', 'chr2:220285666', 'chr3:14174427',
              'chr1:201334382', 'chr10:75865065', 'chr12:22068849', 'chr5:155771579', 'chr3:38592406', 'chr6:7585967']

    conn = sqlite3.connect(db)
    c = conn.cursor()
    
    samples_all = c.execute("SELECT SAMPLE, SERIE FROM snpcheck").fetchall()
    
    sql = f"""SELECT DATA FROM snpcheck 
              WHERE (SAMPLE='{sample_to_check}' 
              AND SERIE ='{serie_to_check}')
              """
    
    snpchecks_to_check = json.loads(c.execute(sql).fetchone()[0])
    snpcheck_ngs_to_check = snpchecks_to_check['NGS']

    for sample, serie in samples_all:
        snpchecks = json.loads(c.execute(f"SELECT DATA FROM snpcheck WHERE (SAMPLE='{sample}' AND SERIE='{serie}')").fetchone()[0])
        snpcheck_alt = snpchecks['ALT']
        snpcheck_alt_calls = list()
        snpcheck_ngs_calls = list()
        
        for locus in alt_loci:
            if snpcheck_alt[locus] in ['WT', 'HET', 'HOM']:
                snpcheck_alt_calls.append(snpcheck_alt[locus])
                snpcheck_ngs_calls.append(snpcheck_ngs_to_check[locus])
        if snpcheck_alt_calls == snpcheck_ngs_calls and len(snpcheck_alt_calls) > 9:
            print(sample_to_check, sample, serie)
            

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Delete sample from database")
    parser.add_argument('-db', '--database', type=str, required=True,
                        help='BED-file to check')
    parser.add_argument('-sa', '--sample', type=str, required=True,
                        help='Sample to delete')
    parser.add_argument('-se', '--serie', type=str, required=True,
                        help='Serie for sample')

    args = parser.parse_args()
    database = args.database
    sample = args.sample
    serie = args.serie
    find_snpcheck(sample, serie, database)

