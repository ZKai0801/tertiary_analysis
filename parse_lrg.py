__doc__ = """
parse LRG database in advance, to facilitate downstream analysis (see annotate_lrg.py)
output file is a tab-delimited file, arranging like this:
    gene_name   transcript_name   LRG_record   strand
    ...

Usage:
[kai@admin]$ python3 parse_LRG_database.py [LRG_GRCh37.bed] [output_filename]

The LRG database was downloaded from https://www.lrg-sequence.org/data/

update:
    Also parse refFlat file to get transcript with longest length
"""
__author__ = "Kai"
__date__ = "2019/10/31"
__version__ = 0.2


import sys
import os
import re
import argparse
import pandas as pd


def parse_lrg(lrg, ofname):
    LRG_genes = {}
    LRG_transcripts = {}
    with open(lrg, "r") as fh:
        for line in fh:
            if line.startswith('track'):
                continue

            match = re.search(r"(LRG_\d+)\((\S+)\)", line.split()[3])
            if match:
                LRG_genes[match.group(2)] = (match.group(1), line.split()[5])
                continue

            match = re.search(r'(LRG_\d+)t1\((N[MR]_\S+)\)', line.split()[3])
            if match:
                LRG_transcripts[match.group(1)] = match.group(2).split("|")[0]
            
           
    # print some logging info 
    print('** The number of genes identified in LRG: {}'.format(len(LRG_genes)))
    print('** The number of transcripts identified in LRG: {}'.format(len(LRG_transcripts)))
    
    # write to output file
    with open(ofname, "w") as fh:
        for gene in LRG_genes:
            fh.write("{}\t{}\t{}\n".format(gene, LRG_transcripts[LRG_genes[gene][0]], LRG_genes[gene][0]))
    print("** Output has been written to {}".format(ofname))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("LRG", help="the LRG database")
    parser.add_argument("-o", "--output", default="clinic_transcript.tsv", help="the output filename")
    args = parser.parse_args()

    parse_lrg(args.LRG, args.output)
    
