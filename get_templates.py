#!/usr/bin/env python3

import requests
from pathlib import Path

def make_output_folder():

    '''make folder to store ouput'''

    output_path = Path('outputs')

    if not output_path.exists():
        output_path.mkdir()


def blast_seq(fasta):

    '''blast a target sequence to identify homologs in the PDB

    Args:
        fasta (str) : file with amino acid sequence for target protein

    Returns:
        blast_txt (str) : blast results text
    '''

    # output dest


    # blast and save results


    # print hits
    pass



def download_sequences(pdb_list):

    pass


def clean_structure(pdb, chain):

    pass


def download_structures(pdb_list):

    pass

def remove_dups(fasta):

    # drop by uniprot id, since raw sequence will not be clean
    pass

def seq_align_hits(fasta, n_templates):

    pass

def write_grishin(fasta):

    pass

def write_ros_xml(grishin):

    pass

def write_flags():

    pass

def write_submission():

    pass
