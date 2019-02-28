#!/usr/bin/env python
import os
import json
import argparse
import subprocess

from ngsscriptlibrary import TargetDatabase

HOME = os.path.expanduser('~')
REPO = os.path.join(HOME, 'Desktop', 'ngstargets')
DB = os.path.join(REPO, 'varia', 'captures.sqlite')

def check_if_file_exists(fn):
    if not os.path.isfile(fn):
        raise ValueError(f'{fn} does not exist')

def copy_file_to_repo(file_location, target_type, repo=None):
    if repo is None:
        repo = REPO
    if target_type == 'capture':
        command = f'cp {file_location} {repo}/captures'
    elif target_type == 'pakket':
        command = f'cp {file_location} {repo}/pakketten'
    elif target_type == 'panel':
        command = f'cp {file_location} {repo}/panels'
    else:
        raise ValueError('Target type moet capture, pakket of panel zijn')
    subprocess.call(command, shell=True, stderr=subprocess.STDOUT)
    return


def get_genelist(genefile):
    with open(genefile) as f:
        genelist = json.load(f)
        genelist = json.dumps(genelist)
        return genelist


def get_info_for_genesis(genesis):
    TD = TargetDatabase(DB)
    info = TD.get_info_for_genesis(genesis)
    return info


def get_size_targts(bed):
    size_total = 0
    with open(bed) as f:
        for line in f:
            chromosome, start, end, *_ = line.split()
            size_target = int(end) - int(start)
            size_total += size_target
    return size_total


def get_new_capture(current_capture):
    capture_name, current_capture_version = current_capture.split('v')
    new_capture_version = int(current_capture_version) + 1
    new_capture = f'{capture_name}v{new_capture_version}'
    return new_capture
    
def add_capture(genesis_info, location, oid):    
    if not oid:
        parser.error('OID required for --capture')

    current_capture = genesis_info['capture']
    new_capture = get_new_capture(current_capture)

    genefile = f'{location}/{new_capture}_genes.json'
    picard_intervals = f'{location}/{new_capture}_target.interval_list'
    capture_bed = f'{location}/{new_capture}_target.bed'
    annotated_bed = f'{location}/{new_capture}_target.annotated'
    generegion_bed = f'{location}/{new_capture}_generegions.bed'

    files_to_copy = [picard_intervals, capture_bed,
                     annotated_bed, generegion_bed]

    check_if_file_exists(genefile)

    for fn in files_to_copy:
        check_if_file_exists(fn)

    for fn in files_to_copy:
        copy_file_to_repo(fn, 'capture')

    genes = get_genelist(genefile)
    size = get_size_targts(capture_bed)
    print(new_capture, size, genes)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add BED and info to repo")
    parser.add_argument('-g', '--genesis', type=str, required=True,
                        help='Genesis code')
    parser.add_argument('-l', '--location', type=str, required=True,
                        help='Folder with required files')
    parser.add_argument('--oid', type=str,
                        help='Capture OID')
    parser.add_argument('--capture', action='store_true',
                        help='Add capture to repo')
    parser.add_argument('--pakket', action='store_true',
                        help='Add pakket to repo')
    parser.add_argument('--panel', action='store_true',
                        help='Add panel to repo')

    args = parser.parse_args()

    if not args.capture and not args.pakket and not args.panel:
        parser.error('--capture, --pakket of --panel is required.')

    genesis = args.genesis
    location = args.location
    new_capture = args.capture
    new_pakket = args.pakket
    new_panel = args.panel
    oid = args.oid

    genesis_info = get_info_for_genesis(genesis)

    if new_capture:
        add_capture(genesis_info, location, oid)
