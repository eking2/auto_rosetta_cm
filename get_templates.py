#!/usr/bin/env python3

from Bio import SeqIO, SearchIO
import requests
import pandas as pd
from pathlib import Path
from io import StringIO
from bs4 import BeautifulSoup
import subprocess
import os
import shutil
from jinja2 import Environment, FileSystemLoader

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

    # check if blast df already written
    name = blast_xml.split('_')[0].split('/')[-1]
    blast_df_path = f'outputs/{name}_blast_df.csv'
    if Path(blast_df_path).exists():
        print(f'{name} already to blast df')
        return

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

    df['uniprot'] = df['pdb_id'].apply(get_uniprot)
    df.to_csv(blast_df_path, index=False)


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
    # no space between pdb and chain because of partial align requirement for 5 chars
    out_path = f'outputs/{pdb}{chain}.pdb'

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

    '''select templates from blast search, download amino acid sequences, and mafft align

    Args:
        fasta (str) : file with template amino acid sequence
        blast_df (str) : blast dataframe csv
        n_templates (int) : number of templates to use

    Returns:
        templates (list) : filenames for template pdbs
    '''

    # check if aligned fasta already written
    name = fasta.split('.')[0]
    hits_out_path = f'outputs/{name}_hits.fasta'
    aln_out_path = f'outputs/{name}_aligned.fasta'


    # select n_templates to download, drop duplicates first
    df = pd.read_csv(blast_df)

    # keep highest seq id and query cov
    df = df.sort_values(by=['seq_ident', 'query_cov'], ascending=False)
    df = df.drop_duplicates(subset='uniprot')
    df = df.iloc[:n_templates]

    # output selected
    templates = []

    for idx, row in df.iterrows():
        templates.append(f"{row['pdb_id'].upper()}{row['chains'][0]}.pdb")

    if Path(aln_out_path).exists():
        print(f'already aligned {name}')
        return templates

    # download structures and sequences
    records = []

    # append template
    records.append(SeqIO.read(fasta, 'fasta'))

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

    return templates


# can combine with seq_align_hits
def seq_align_templates(fasta, templates):

    '''download amino acid sequences for input templates and mafft align

    Args:
        fasta (str) : file with template amino acid sequence
        templates (list) : names of templates and chain (ie. 4E5NA.pdb)
    '''

    # check if aligned fasta already written
    name = fasta.split('.')[0]
    hits_out_path = f'outputs/{name}_hits.fasta'
    aln_out_path = f'outputs/{name}_aligned.fasta'

    if Path(aln_out_path).exists():
        print(f'already aligned {name}')
        return

    # download structures and sequences
    records = []

    # append template
    records.append(SeqIO.read(fasta, 'fasta'))

    for template in templates:

        pdb = template.split('_')[0]
        chain = template.split('_')[1]

        record = download_sequence(pdb)
        record.id = f'{pdb}_{chain}'
        record.description = ''

        clean_structure(pdb, chain)
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

    '''convert msa to grishin format for rosetta partial threading

    Args:
        fasta (str) : file with msa of hits
    '''

    grish_template = '''
## {target} {template}_thread
#
scores_from_program: 0
0 {target_seq}
0 {template_seq}
--
    '''

    # get name of target and sequence
    records = list(SeqIO.parse(fasta, 'fasta'))
    target = records[0].name.replace('_', '')
    target_seq = records[0].seq

    # loop over templates, align to target
    # skip first, target
    out = []
    for i in range(1, len(records)):
        template = records[i].name.replace('_', '')
        template_seq = records[i].seq
        grish = grish_template.format(target = target, template=template, target_seq = target_seq, template_seq=template_seq)
        out.append(grish)

    grish_out = ''.join(out)
    name = fasta.split('_aligned')[0].split('/')[-1]
    Path(f'outputs/{name}.gri').write_text(grish_out)


def run_partial_thread(fasta, grish_aln, templates):

    '''run partial threading of aligned sequence onto template structures

    Args:
        fasta (str) : path for target protein fasta
        grish_aln (str) : grishin format alignment for target to templates
        templates (list) : list of pdb files for templates
    '''

    # list to str
    templates = [f'outputs/{x.replace("_", "")}' for x in templates]
    templates = ' '.join(templates)

    flags = f'''-in:file:fasta {fasta}
-in:file:alignment {grish_aln}
-in:file:template_pdb {templates}
-ignore_unrecognized_res 1'''

    flags_path = f'outputs/partial_thread_flags'
    Path(flags_path).write_text(flags)

    # run partial thread
    # must export ROSETTA3 to install dir before
    env = os.environ.copy()
    ros = env['ROSETTA3']
    subprocess.run([f'{ros}/bin/partial_thread.static.linuxgccrelease', f'@{flags_path}'])

    # move threaded pdbs to outputs
    threaded = Path().glob('*thread.pdb')
    for pdb in threaded:
        shutil.move(str(pdb), 'outputs')


def write_hybridize_flags(fasta, nstruct=5):

    '''write hybridize flag file and shell script to run rosetta script

    Args:
        fasta (str) : path for target protein fasta
        nstruct (int) : number of modeling trajectories
    '''

    flags = f'''-default_max_cycles 200
-nstruct {nstruct}
-dualspace
-in:file:fasta {fasta}
-parser:protocol hyb_run.xml'''

    name = Path(fasta).stem
    Path(f'outputs/{name}_hyb_flags').write_text(flags)

    shell_cmd = f'$ROSETTA3/bin/rosetta_scripts.static.linuxgccrelease @{name}_hyb_flags'
    Path('outputs/run_hybridize.sh').write_text(shell_cmd)


def write_hybridize_xml(templates):

    '''write rosetta script xml input with templates pdb files

    Args:
        templates (list) : list of template pdb files
    '''

    # add thread suffix
    for i, template in enumerate(templates):
        name = template.split('.')[0]
        name = name.replace('_', '')
        templates[i] = f'{name}_thread.pdb'

    # replace with jinja
    env = Environment(loader=FileSystemLoader('input_templates'))
    template = env.get_template('hyb.xml')

    output = template.render(templates=templates)
    Path('outputs', 'hyb_run.xml').write_text(output)

