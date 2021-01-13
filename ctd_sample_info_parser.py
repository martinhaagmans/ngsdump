import csv
import sqlite3

DB = '/data/dnadiag/databases/patientinfo.sqlite'

def parse_input(input_file):
    with open(input_file) as f:
        reader = csv.reader(f)
        _header = next(reader)
        cnv_info = list()
        patient_info = list()
        for line in reader:
            serie, sample, material, sex, request, dob, deadline = line
            patient_info.append([serie, sample, sex, request, dob, deadline])
            if not material.lower() == 'bld':
                cnv_info.append([serie, sample, material])
    return (cnv_info, patient_info)


def create_cnv_info_file(cnv_info, file_dir):
    cnv_info_file = os.path.join(file_dir, 'cnv_info.txt')
    with open(cnv_info_file, 'w') as f:
        for serie, sample, material in cnv_info:
            f.write('{}\t{}\t{}\n'.format(serie, sample, material))
    return cnv_info_file


def add_patient_info_to_db(patient_info):
    conn = sqlite3.connect(DB)
    c = conn.cursor()
    for serie, sample, sex, request, dob, deadline in patient_info:
        sql = """INSERT INTO patientinfo 
        (SERIE, SAMPLE, SEX, FF, DOB, DEADLINE)
        VALUES ("{}", "{}", "{}", "{}", "{}", "{}")
        """.format(serie, sample, sex, request, dob, deadline)
        c.execute(sql)
        conn.commit()
    conn.close()
    return


if __name__ == '__main__':
    import os
    import argparse
    from ngsscriptlibrary import run_command

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="Input file", required=True)

    args = parser.parse_args()
    cnv_info, patient_info = parse_input(args.input)
    
    add_patient_info_to_db(patient_info)

    if cnv_info:
        file_dir = os.path.dirname(args.input)
        cnv_info_file = create_cnv_info_file(cnv_info, file_dir)
        run_command('CNV --create -c BLA -i {}'.format(cnv_info_file))
    
