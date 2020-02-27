import os 
import json
import sqlite3
import pandas as pd

def get_archive(db, capture):
    """Get coveragedata from table and return df"""

    sql = f"SELECT * FROM {capture}"
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute(sql)
    archive = dict()
    for sample, serie, data in c.fetchall():
        archive[sample] = dict()

        archive[sample]['data'] = dict()
        archive[sample]['serie'] = serie

        for tup in json.loads(data):
            target, coverage = tup
            archive[sample]['data'][target] = coverage

    conn.close()
    return archive
    
def archive_to_table(archive):    
    dflist = list()

    for sample in archive.keys():
        dftmp = pd.Series(archive[sample]['data'], name=sample, dtype=float)
        dftmp['serie'] = archive[sample]['serie']
        dftmp['sample'] = sample
        dftmp = dftmp.transpose()
        dflist.append(dftmp)

    df = pd.concat(dflist, axis=1)
    df = df.transpose().set_index('sample').sort_index().sort_index(axis=1)
    return df
    
def get_annot_table(db, capture):
    annot = pd.read_sql(f'SELECT * FROM {capture}annot', con = sqlite3.connect(db))
    annot = annot.set_index('target')
    return annot

def main(db, capture):
    archive = get_archive(db, capture)
    df =  archive_to_table(archive)
    dfa = get_annot_table(db, capture)
    df = df.transpose().join(dfa)
    cols = list(df.columns)
    cols = cols[-1:] + cols[:-1]
    ix = list(df.index)
    ix = ix[-1:] + ix[:-1]
    df = df.reindex(ix) 
    df = df.reindex(cols, axis=1)
    df.index = [_.replace('_', ':', 1) for _ in df.index]
    df.index = [_.replace('_', '-', 1) for _ in df.index]
    df.to_excel(f'{capture}.xlsx')
    return
    
if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-db", "--database", type=str, 
                        help="Database", required=True)
    parser.add_argument("-c", "--capture", type=str, 
                        help="Database", required=True)                        
    args = parser.parse_args()
    
    db = args.database
    capture = args.capture
    main(db, capture)
    
