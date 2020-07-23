#!/usr/bin/env python3

import argparse
import get_templates

def parse_args():

    parser = argparse.ArgumentParser()

    # read template fasta
    parser.add_argument('-f', '--fasta', type=str, required=True,
                        help='Fasta file with template protein amino acid sequence')

    # run blast or setup cm
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-b', '--blast', action='store_true',
                        help='Run blast to identify homologs')
    group.add_argument('-c', '--cm', action='store_true',
                        help='Setup RosettaCM for homology modeling')

    # cm options
    parser.add_argument('-n', '--num_templates', type=int, default=5,
                        help='Number of templates to use for homology modeling')
    parser.add_argument('-l', '--low_id', type=float, default=0.2,
                        help='Sequence identity cutoff, drop hits below')
    parser.add_argument('-q', '--query_over', type=float, default=0.5,
                        help='Minimum query overlap, drop hits below')
    parser.add_argument('-t', '--templates', nargs='*', default=[],
                        help='Target PDB templates to use, format as XXXX_Chain (ie. 4e5n_A)')

    return parser.parse_args()
    

if __name__ == '__main__':

    args = parse_args()
    print(args)

    # make output folder
    get_templates.make_output_folder()

    # setup rosetta cm
    if args.cm:

        # use target templates
        if len(args.templates) > 0:
            pass

        pass

    # run blast search on rcsb
    else:
        pass


