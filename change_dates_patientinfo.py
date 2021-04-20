import sqlite3

conn = sqlite3.connect('patientinfo.sqlite')
c = conn.cursor()
c.execute('SELECT * FROM patientinfo')
all_info = c.fetchall()
conn.close()

out = list()
for info in all_info:
    serie, dnr, sex, ff, dob = info
    if serie > 745 and serie < 756:
        d, m, y = dob.split('-')
        if int(d) <= 12 and dob != '0-1-1900' :
            dob = f'{m}-{d}-{y}'
    out.append((serie, dnr, sex, ff, dob))

conn = sqlite3.connect('patientinfo2.sqlite')
c = conn.cursor()
sql =""" CREATE TABLE 'patientinfo' ( 
         'SERIE' INTEGER NOT NULL, 
         'SAMPLE' TEXT NOT NULL, 
         'SEX' TEXT, 
         'FF' TEXT, 
         'DOB' TEXT, 
         PRIMARY KEY('SERIE','SAMPLE') 
         )"""
         
c.execute(sql)

for serie, dnr, sex, ff, dob in out:
    c.execute(f"INSERT INTO patientinfo VALUES ({serie}, '{dnr}', '{sex}', '{ff}', '{dob}')")

conn.commit()
conn.close()


