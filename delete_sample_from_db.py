import argparse
import sqlite3


class DeleteSample:
    def __init__(self, database, sample, serie):
        self.sample = sample
        self.serie = serie
        self.conn = sqlite3.connect(database)
        self.c = self.conn.cursor()

    def __del__(self):
        self.conn.close()

    def _get_tables(self):
        self.c.execute("SELECT name FROM sqlite_master WHERE type='table'")
        tables = [val for tup in self.c.fetchall() for val in tup]
        return tables
    
    def delete(self):
        for table in self._get_tables():
            sql = f"""DELETE FROM {table}
                WHERE (SAMPLE='{self.sample}' 
                AND SERIE='{self.serie}')
                """
            self.c.execute(sql)
            self.conn.commit()


    


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

    DeleteSample(database, sample, serie).delete()
