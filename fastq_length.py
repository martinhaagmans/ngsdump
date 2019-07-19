import argparse


def get_readlengths(fastq_file):
    read_lengths = list()
    
    prefix = fastq_file.split('.')[0]
    outfile = f'{prefix}_lenghts.txt'
    
    with open(fastq_file) as f_in, open(outfile, 'w') as f_out:
        seq = str()
        for line in f_in:
            line = line.rstrip('\n')
            if line.startswith('@'):
                if seq:
                    f_out.write(f'{name}\t{str(len(seq))}\n')
                    seq = str()
                name = line.split(' ')[0]
            else:
                seq = line
                read_lengths.append(len(seq)) 
        f_out.write(f'{name}\t{str(len(seq))}\n')
    return read_lengths

   
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq", type=str, 
                        help="Fastq file", required=True)
    args = parser.parse_args()                          
    get_readlengths(args.fastq)

