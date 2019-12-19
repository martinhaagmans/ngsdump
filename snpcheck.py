import os
import sys
import json
import logging
import sqlite3


class SNPcheck:

    def __init__(self, sample_to_check, db):
       if not os.path.isfile(db):
        logging.error(f'{db} does not exist')
        sys.exit()
       self.conn = sqlite3.connect(db)
       self.c = self.conn.cursor()
       self.sample_to_check = sample_to_check
       self.rs_gpos = { 
        'rs1071646': 'chr15:63351840',
        'rs4485000': 'chr10:18789724',
        'rs1063243': 'chr7:91726927',
        'rs2250736': 'chr3:53700550',
        'rs12920': 'chr2:220285666',
        'rs4685076': 'chr3:14174427',
        'rs3729547': 'chr1:201334382',
        'rs767809': 'chr10:75865065',
        'rs3759236': 'chr12:22068849',
        'rs1801193': 'chr5:155771579',
        'rs1805126': 'chr3:38592406',
        'rs2744380': 'chr6:7585967',
        'chr15:63351840': 'rs1071646',
        'chr10:18789724': 'rs4485000',
        'chr7:91726927': 'rs1063243',
        'chr3:53700550': 'rs2250736',
        'chr2:220285666': 'rs12920',
        'chr3:14174427': 'rs4685076',
        'chr1:201334382': 'rs3729547',
        'chr10:75865065': 'rs767809',
        'chr12:22068849': 'rs3759236',
        'chr5:155771579': 'rs1801193',
        'chr3:38592406': 'rs1805126',
        'chr6:7585967': 'rs2744380'
        }
       

    def __del__(self):
        try:
            self.conn.close()
        except AttributeError:
            pass
        
    def get_all_samples(self):
        self.c.execute('SELECT SAMPLE FROM snpcheck')
        samples =  [val for tup in self.c.fetchall() for val in tup]
        return samples

    def get_all_snpchecks(self):
        samples = self.get_all_samples()
        all_snpchecks = dict()
        for sample in samples:
            self.c.execute(f'SELECT DATA FROM snpcheck WHERE SAMPLE="{sample}"')
            snpcheck = self.c.fetchone()[0]
            snpcheck = json.loads(snpcheck)  
            all_snpchecks[sample] = snpcheck
        return all_snpchecks
        
    def get_all_ngs_snpcheck_strings(self):
        snpcheck = self.get_all_snpchecks()
        samples = snpcheck.keys()
        all_ngs_snpcheck_strings = dict()
        for sample in samples:
            ngs_snpcheck_string = list()
            for _locus, call in sorted(snpcheck[sample]['NGS'].items()):
                ngs_snpcheck_string.append(call)
            ngs_snpcheck_string = ''.join(ngs_snpcheck_string)
            all_ngs_snpcheck_strings[sample] = ngs_snpcheck_string
        return all_ngs_snpcheck_strings
              
    def get_all_alt_snpcheck_strings(self):
        snpcheck = self.get_all_snpchecks()
        samples = snpcheck.keys()
        all_alt_snpcheck_strings = dict()
        for sample in samples:
            alt_snpcheck_string = list()
            for _locus, call in sorted(snpcheck[sample]['ALT'].items()):
                alt_snpcheck_string.append(call)
            alt_snpcheck_string = ''.join(alt_snpcheck_string)
            all_alt_snpcheck_strings[sample] = alt_snpcheck_string
        return all_alt_snpcheck_strings
          
    def compare_ngs_snpchecks(self):
        has_same_ngs_snpcheck = set()
        all_snpchecks = self.get_all_ngs_snpcheck_strings()
        try:
            snpcheck_to_check = all_snpchecks[self.sample_to_check]
        except KeyError:
            logging.error(f'{self.sample_to_check} not in db')
            sys.exit()
        for sample, snpcheck in all_snpchecks.items():
            if snpcheck_to_check == snpcheck:
                has_same_ngs_snpcheck.add(sample)
                
        return has_same_ngs_snpcheck

   


if __name__ == '__main__':
    import argparse
    
    logging.basicConfig(format='%(levelname)s: %(message)s')
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample", type=str, 
                        help="SampleID", required=True)
    parser.add_argument("-db", "--database", type=str, 
                        help="Metrics database", required=True)                        
                        
    args = parser.parse_args()
    print(SNPcheck(args.sample, args.database).compare_ngs_snpchecks())
