import os
import sqlite3

DB = '/data/dnadiag/databases/patientinfo.sqlite'


def read_patientinfo(patientinfo):
    data = list()
    _folder, fn = os.path.split(patientinfo)
    serie = fn.split('_')[0].replace('MS', '')
    with open(patientinfo) as f:
        for line in f:
            if line not in ['\n', '\r\n']:
                data.append('{}\t{}'.format(serie, line))
    return data


def add_to_db(data):
    conn = sqlite3.connect(DB)
    c = conn.cursor()    
    for line in data:
        serie, sample, sex, ff, dob = line.split()
        try:
            c.execute('INSERT INTO patientinfo VALUES ("{}", "{}", "{}", "{}", "{}")'.format(serie, sample, sex, ff, dob))
        except sqlite3.IntegrityError as e:
            print(e)
        else:
            conn.commit()
    conn.close()
    return


def main(patientinfo):
    data = read_patientinfo(patientinfo)
    add_to_db(data)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, 
                        help="File met patient info", required=True)
    args = parser.parse_args()
    main(args.input)