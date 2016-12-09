import pandas as pd
import csv
import datetime

excelfile = 'MiSeqExp161(0-lijsttotMiSeqrun)_20161031.xlsx'
serie = 'MSTEST'
readlength = 151
now = datetime.datetime.now()
date = now.strftime("%m/%d/%Y")

#adapter files
adapterfile = 'adapters.txt'
adaptercodefile = 'adaptercodes.txt'

def create_adapter_dict(adapterfile):
    """Read tab seperated file with: number adapterseq1 adapterseq2.
    Return dict with number as key and tuple of sequences as value.
    """
    adapterdict = dict()
    with open(adapterfile) as f:
        for line in f:
            if not line.split():
                continue
            else:
                index, left, right = line.strip().split('\t')
                adapterdict[int(index)] = [left, right]
    return adapterdict

def create_reverse_adapter_dict(adapterfile):
    """Read tab seperated file with: number adapterseq1 adapterseq2.
    Return dict with number as key and tuple of sequences as value.
    """
    reverse_adapterdict = dict()
    with open(adapterfile) as f:
        for line in f:
            if not line.split():
                continue
            else:
                index, left, right = line.strip().split('\t')
                reverse_adapterdict[left, right] = [int(index)]
    return reverse_adapterdict

def create_adapter_code_dict(adaptercodefile):
    """Read tab seperated file with adapterseq adapterID.
    Return dict with sequence as key and ID as value.
    """
    adapterdict = dict()
    with open(adaptercodefile) as f:
        for line in f:
            if not line.split():
                continue
            else:
                code, sequence = line.strip().split('\t')
                adapterdict[sequence] = code 
    return adapterdict

def parse_excelfile(excelfile):
    """Read excel file. Return dataframe with 2 columns."""
    df = pd.read_excel(excelfile, header=5)
    df = df[['Dnummer', 'Gen']]
    df.columns = ['Dnummer', 'Capture']
    return df

def pick_barcodes(df):
    """Read dataframe.
    Return dict with number as key and and sequence
    """
    patients = df['Dnummer'].values
    barcode = 97
    patd = dict()
    adaptersequences = create_adapter_dict(adapterfile)

    for i in range(len(patients)):
        if barcode < 109:
            first_barcode = barcode
        elif barcode > 192:
            barcode = first_barcode + 1
        patd[patients[i]] = adaptersequences[barcode]
        barcode += 13
    return patd

def writelist(list_to_print, f):
    """Prints items in list to fileobject.
    If item > 1: print csv tuple
    """
    for val in list_to_print:
        if not isinstance( val, int ):
            if '[' in val and ']' in val:
                f.write(str(val))
            elif len(val) == 2:
                f.write('{} ,'.format(val[0]))
                f.write(str(val[1]))
        elif isinstance( val, int ):
            f.write('{}'.format(val))
        f.write('\n')

def main():
    # Create dicts with adapter info
    
    adaptercodes = create_adapter_code_dict(adaptercodefile)
    reverse_adapterdict = create_reverse_adapter_dict(adapterfile)

    df = parse_excelfile(excelfile)

    barcodes = pick_barcodes(df)

    out = [(k,k, '', '', adaptercodes[v[1]], v[1], adaptercodes[v[0]], v[0], serie)
           for k,v in barcodes.items()]

    dfout = pd.DataFrame(out)
    dfout = dfout.join(df, rsuffix='right')

    dfout.drop('Dnummer', axis=1, inplace=True)

    dfout.columns = ['Sample_ID', 'Sample_Name', 'Sample_Plate', 'Sample_Well',
                     'I7_Index_ID','index', 'I5_Index_ID', 'index2',
                     'Sample_Project', 'Description']

    captures = '_'.join(dfout['Description'].unique())


    #Define samplesheet
    header = ['[Header]',
              ['IEMFileVersion', 4] ,
              ['Investigator Name', 'Robert'],
              ['Experiment Name', serie],
              ['Date', date],
              ['Workflow', 'GenerateFASTQ'],
              ['Application', 'FASTQ Only'],
              ['Assay', 'BIOO set 2'],
              ['Description', captures],
              ['Chemistry', 'Amplicon']]


    reads = ['[Reads]', readlength, readlength]

    settings = ['[Settings]',
                ['ReverseComplement', 0],
                ['Adapter', 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'],
                ['AdapterRead2', 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT']]

    samplesheet = 'samplesheet.csv'


    #Output
    with open('barcode_pipetteerschema.txt', 'w') as f:
        [[f.write('{}\t{}\n'.format(k, i)) for i in reverse_adapterdict[(v[0], v[1])] ]
          for k,v in barcodes.items()]


    with open(samplesheet, 'w') as f:
        writelist(header, f)
        f.write('\n')
        writelist(reads, f)
        f.write('\n')
        writelist(settings, f)
        f.write('\n')
        f.write('[Data]\n')
        dfout.to_csv(f, index=None)
        
        
if __name__ == "__main__":
    main()


























