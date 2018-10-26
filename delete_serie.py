class DeleteSerie:
    def __init__(self, db):
        self.db = db
        self.conn = sqlite3.connect(db)
        self.c = self.conn.cursor()

    def __del__(self):
        self.conn.close()
        

    def get_table_names(self):
        self.c.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [val for tup in self.c.fetchall() for val in tup]
        return tables
        
    def check_if_serie_exists(self, table, serie):
        self.c.execute(f"SELECT count(*) FROM {table} WHERE serie = '{serie}'")
        rows_for_serie = self.c.fetchone()[0]
        return rows_for_serie != 0

    def delete_serie_from_table(self, table, serie):
        if not self.check_if_serie_exists(table, serie):
            raise ValueError(f'{serie} niet in {table}')        
        else:
            self.c.execute(f"DELETE FROM {table} WHERE SERIE='{serie}'")
            self.conn.commit()
            return

    def delete_serie_from_db(self, serie):
        for table in self.get_table_names():
            try:
                self.delete_serie_from_table(table, serie)
            except ValueError as e:
                print(e)
            else:
                print(f"Deleted {serie} from {table} in {self.db}")
    

x = DeleteSerie('Desktop/metrics.sqlite')
x.delete_serie_from_db(444)



