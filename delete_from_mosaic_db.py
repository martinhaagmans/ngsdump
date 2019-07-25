import argparse
import sqlite3

def delete_sample_from_mosaic_database(sample, db, table=None):
    if table is None:
        table = 'patients'

    sql = f'DELETE FROM {table} WHERE SAMPLE="{sample}"'

    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute(sql)
    conn.commit()
    c.execute('vacuum')
    conn.close()
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--sample", type=str, 
                        help="Sample to remove", required=True)
                        
    parser.add_argument("-db", "--database", type=str, 
                        help="Database to remove sample from", required=True)
                        
                                               
    args = parser.parse_args()                          
    sample = args.sample
    db = args.database

    delete_sample_from_mosaic_database(sample, db, table=None)

