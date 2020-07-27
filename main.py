#!/usr/bin/env python3

import argparse
from pathlib import Path
from get_templates import *

def parse_args():

    parser = argparse.ArgumentParser()

    # read template fasta
    parser.add_argument('-f', '--fasta', type=str, required=True,
                        help='Fasta file with template protein amino acid sequence')

    # separate args to run blast search and cm
    parser.add_argument('-b', '--blast', action='store_true',
                        help='Run blast to identify homologs')

    # cm options
    subparsers = parser.add_subparsers(dest='process')
    cm_parser = subparsers.add_parser('cm', help='Setup RosettaCM for homology modeling')

    cm_parser.add_argument('-n', '--num_templates', type=int, default=5,
            help='Number of templates to use for homology modeling (default: 5)')
    cm_parser.add_argument('-d', '--decoys', type=int, default=500,
            help='Number of decoys to generate (default: 500)')
    cm_parser.add_argument('-t', '--templates', nargs='*', default=[],
                        help='Target PDB templates to use, format as XXXX_Chain (ie. 4e5n_A)')

    return parser.parse_args()


def run_blast(fasta):

    '''run blast search'''

    name = Path(fasta).stem
    blast_seq(args.fasta)
    blast_xml = f'outputs/{name}_blast.xml'
    parse_blast_xml(blast_xml)

def run_cm_setup(fasta, templates, n_templates, n_decoys):

    '''perform blast search, select templates, and write inputs for rosetta cm'''

    # auto template selection
    if len(templates) == 0:

        # start with blast to select templates
        name = Path(fasta).stem
        run_blast(fasta)
        blast_df = f'outputs/{name}_blast_df.csv'

        # msa hits
        templates = seq_align_hits(fasta, blast_df, n_templates)

    # manual template selection
    else:
        # download seq+pdb and align sequences
        seq_align_templates(fasta, templates)

        # add file suffix if missing, for writing grishin
        for i, template in enumerate(templates):
            if not template.endswith('pdb'):
                templates[i] = f'{template}.pdb'


    # partial thread and write inputs
    aln_fasta = f'outputs/{name}_aligned.fasta'
    grish_aln = f'outputs/{name}.gri'

    print(f"templates: {', '.join(templates)}")
    write_grishin(aln_fasta)
    run_partial_thread(fasta, grish_aln, templates)
    write_hybridize_flags(fasta, n_decoys)
    write_hybridize_xml(templates)


if __name__ == '__main__':

    args = parse_args()

    # make output folder
    make_output_folder()

    # blast search
    if args.blast:
        print(f'running blast search on: {args.fasta}')
        run_blast(args.fasta)

    # setup rosetta cm
    if args.process == 'cm':
        print(f'setting up rosetta cm: {args.fasta}')
        run_cm_setup(args.fasta, args.templates, args.num_templates, args.decoys)
