#!/usr/bin/env python3

from Bio import SeqIO, SearchIO
import requests
import pandas as pd
from pathlib import Path
from io import StringIO
from bs4 import BeautifulSoup
import subprocess

'''
run blast search to identify homologs

1. read template fasta
2. run blast through rcsb rest api

prepare input files for RosettaCM homology modeling

1. download template pdbs and fastas
2. clean pdbs
3. align sequences
4. make grishin
5. partial thread
6. write rosetta xml input 
7. hybridize

'''

def make_output_folder():

    '''make folder to store ouput'''

    output_path = Path('outputs')

    if not output_path.exists():
        output_path.mkdir()


def blast_seq(fasta):

    '''blast a target sequence to identify homologs in the PDB.

    Args:
        fasta (str) : file with amino acid sequence for target protein
    '''
    
    name = Path(fasta).stem
    out_name = f'outputs/{name}_blast.xml'

    # do not blast if output already exists
    if Path(out_name).exists():
        print(f'{name} already blasted')
        return

    # load fasta
    record = SeqIO.read(fasta, 'fasta')
    seq = record.seq

    # blast
    url = f'https://www.rcsb.org/pdb/rest/getBlastPDB1?sequence={seq}&eCutOff=10.0&matrix=BLOSUM62&outputFormat=XML'
    r = requests.get(url)
    assert r.status_code == 200, 'invalid blast search'

    # save results
    content = r.text
    Path(out_name).write_text(content)


def get_uniprot(pdb):

    '''get uniprot accession from pdb

    Args:
        pdb (str) : 4-letter pdb code

    Returns:
        uniprot_id (str) : uniprot accession id
    '''

    url = f'https://www.rcsb.org/pdb/rest/describeMol?structureId={pdb}'
    r = requests.get(url)
    assert r.status_code == 200, f'invalid pdb: {pdb}'

    soup = BeautifulSoup(r.content, 'lxml')

    try:
        uniprot = soup.find('macromolecule').findChildren('accession')[0]['id']
    except:
        uniprot = None

    return uniprot


def parse_blast_xml(blast_xml):

    '''parse blast xml to dataframe.

    Args:
        blast_xml (str) : blast xml file
    '''

    record = SearchIO.read(blast_xml, 'blast-xml')

    # initial query len
    query_len = record.seq_len

    to_save = []
    for i in range(len(record)):

        query_cov = record.hsps[i].aln_span / query_len
        seq_ident = record.hsps[i].ident_num / query_len

        hit_id = record.hits[i].id
        pdb_id = hit_id[:4]
        chains = hit_id.split('|')[0].split(':')[-1]

        # pdb blast has no desc
        #hit_desc = record.hits[i].description
        # not sure what accession points to
        #hit_accession = record.hits[i].accession

        to_save.append([pdb_id, chains, round(query_cov, 3), round(seq_ident, 3)])

    df = pd.DataFrame(to_save, columns=['pdb_id', 'chains', 'query_cov', 'seq_ident'])

    name = blast_xml.split('_')[0].split('/')[-1]
    df['uniprot'] = df['pdb_id'].apply(get_uniprot)
    df.to_csv(f'outputs/{name}_blast_df.csv', index=False)


def download_sequence(pdb):

    '''download amino acid sequence for target pdb from rcsb.

    Args:
        pdb (str) : 4-character pdb id

    Returns:
        record (BioPython SeqRecord) : fasta record for target pdb
    '''

    url = f'https://www.rcsb.org/fasta/entry/{pdb}'
    r = requests.get(url)
    assert r.status_code == 200, f'invalid request for {pdb}'

    record = SeqIO.read(StringIO(r.text), 'fasta')

    return record


def clean_structure(pdb, chain):

    '''download pdb file and remove waters and ligands, only keep target chain

    Args:
        pdb (str) : 4-character pdb id
        chain (str) : subunit to keep
    '''

    chain = chain.upper()
    out_path = f'outputs/{pdb}_{chain}.pdb'

    if Path(out_path).exists():
        pass

    url = f'https://files.rcsb.org/download/{pdb}.pdb'
    r = requests.get(url)
    assert r.status_code == 200, f'invalid request for pdb {pdb}'

    content = r.text.splitlines()

    to_save = []
    for line in content:
        if line.startswith('ATOM') and line[21] == chain:
            to_save.append(line)

    clean_pdb = '\n'.join(to_save)
    Path(out_path).write_text(clean_pdb)


def seq_align_hits(fasta, blast_df, n_templates):

    '''download amino acid sequences for templates and mafft align

    Args:
        fasta (str) : file with template amino acid sequence
        blast_df (str) : blast dataframe csv
        n_templates (int) : number of templates to use
    '''

    # check if aligned fasta already written
    name = fasta.split('.')[0]
    hits_out_path = f'outputs/{name}_hits.fasta'
    aln_out_path = f'outputs/{name}_aligned.fasta'

    assert not Path(aln_out_path).exists(), f'already aligned {name}'

    # select n_templates to download, drop duplicates first
    df = pd.read_csv(blast_df)

    # keep highest seq id and query cov
    df = df.sort_values(by=['seq_ident', 'query_cov'], ascending=False)
    df = df.drop_duplicates(subset='uniprot')
    df = df.iloc[:n_templates]

    # download structures and sequences
    records = []
    for _, row in df.iterrows():
        chain = row['chains'].split(',')[0]

        # clean record name
        record = download_sequence(row['pdb_id'])
        record.id = f"{row['pdb_id']}_{chain}"
        record.description = ''

        clean_structure(row['pdb_id'], chain)

        # to concat all fasta sequences
        records.append(record)

    # write out hits
    SeqIO.write(records, hits_out_path, 'fasta')

    # mafft align
    proc = subprocess.Popen(['mafft', hits_out_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    out = out.decode('utf8')
    
    # save alignment
    Path(aln_out_path).write_text(out)


def write_grishin(fasta):

    pass

def write_ros_xml(grishin):

    pass

def write_flags():

    pass

def write_submission():

    pass
