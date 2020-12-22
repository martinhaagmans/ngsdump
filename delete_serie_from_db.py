import argparse
import sqlite3


class DeleteSerie:
    def __init__(self, database, serie):
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
                WHERE SERIE='{self.serie}'
                """
            try:
                self.c.execute(sql)
            except sqlite3.OperationalError as e:
                print(table, e)
            else:
                self.conn.commit()
        self.c.execute("VACUUM")


    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Delete serie from database")
    parser.add_argument('-db', '--database', type=str, required=True,
                        help='BED-file to check')
    parser.add_argument('-se', '--serie', type=str, required=True,
                        help='Serie to delete')

    args = parser.parse_args()
    database = args.database
    serie = args.serie

    DeleteSerie(database, serie).delete()
