#!/usr/bin/bash
__doc__ = """
    CIVic is a great database that provides clinical interpretation for variants in cancer:
    https://civicdb.org/home

    This script take VCF file that already been annotated with HGVS nomenclature on each variant
    (contains pHGVS on info section), and output VCF file with annotated civix information 
"""
__author__ = "Kai"
__date__ = "2019/12/30"


import pysam
import argparse
import os
import re
from collections import defaultdict
import pandas as pd
import numpy as np


def main(ifname, ofname, civic):
    civic_variants, civic_phgvs, civic_chgvs, civic_exon = parse_civic(civic)
    with pysam.VariantFile(ifname, "r") as vcfin:
        if "civic" not in vcfin.header.info:
            vcfin.header.info.add("civic", ".", "String", "Information described in CIVic. Format: variant|disease|drugs|evidence_level|evidence_statement|variant_origin|citation_id|citation")
        vcfout = pysam.VariantFile(ofname, "w", header = vcfin.header)
        for rec in vcfin:
            # normalise genomic coordinate 
            if len(rec.ref) == 1 and len(rec.alts[0]) == 1:
                start = rec.pos
                ref, alt = rec.ref, rec.alts[0]

            elif len(rec.ref) > len(rec.alts[0]):
                # deletion
                if len(rec.alts[0]) == 1:
                    start = rec.start + 2
                    ref = rec.ref[1:]
                    alt = "-"
                else:
                    start = rec.start + 2
                    ref = rec.ref[1:]
                    alt = rec.alts[0][1:]
            else:
                # insertion
                if len(rec.alts[0]) == 1:
                    start = rec.start + 1
                    ref = "-"
                    alt = rec.alts[0][1:]
                else:
                    start = rec.start + 1
                    ref = rec.ref[1:]
                    alt = rec.alts[0][1:]
            
            if (rec.chrom.lstrip("chr"), start, ref, alt) in civic_variants:
                rec.info['civic'] = "|".join(civic_variants[(rec.chrom.lstrip("chr"), start, ref, alt)])
                continue

            if rec.info['gene'] in civic_phgvs:
                if rec.info['pHGVS']:
                    phgvs = convert_amino_acid(rec.info['pHGVS'].lstrip("p."))
                    if phgvs in civic_phgvs[rec.info['gene']]:
                        rec.info['civic'] = "|".join(civic_phgvs[rec.info['gene']][rec.info['pHGVS']])
                    if phgvs[:-1] in civic_phgvs[rec.info['gene']]:
                        rec.info['civic'] = "|".join(civic_phgvs[rec.info['gene']][rec.info['pHGVS']])
                    if rec.info['cHGVS'] in civic_chgvs[rec.info['gene']]:
                        rec.info['civic'] = "|".join(civic_chgvs[rec.info['gene']][rec.info['cHGVS']])

            vcfout.write(rec)


def convert_amino_acid(hgvs):
    """
    Convert three-letter abbreviation to one-letter abbreviation, so that hgvs name can meet oncoKB criteria
    e.g. p.Pro373fs --> p.P373fs
    """
    amino_dict = {'Gly': 'G', 'Cys': 'C', 'Glu': 'E', 'Asp': 'D', 'Ile': 'I',
                  'Pro': 'P', 'Tyr': 'Y', 'Lys': 'K', 'Gln': 'Q', 'Trp': 'W',
                  'Leu': 'L', 'Phe': 'F', 'Val': 'V', 'Ser': 'S', 'Met': 'M',
                  'Ala': 'A', 'His': 'H', 'Ter': 'X', 'Asn': 'N', 'Thr': 'T',
                  'Arg': 'R'}
    new_hgvs = hgvs
    for index, letter in enumerate(hgvs):
        if hgvs[index: index+3] not in amino_dict:
            continue
        new_hgvs = new_hgvs.replace(hgvs[index:index+3], amino_dict[hgvs[index:index+3]])
    return new_hgvs


def parse_civic(civic):
    """
    parse civic database into a dictionary. Format:
    civic_database = {gene: {variant: (disease, drug, evidence_level, evidence_statement, variant_origin, citation_id, citation), ...}, ...}
    """
    civic_variants = {}
    civic_phgvs = defaultdict(dict)
    civic_chgvs = defaultdict(dict)
    civic_exon = defaultdict(dict)

    df = pd.read_csv(civic, sep="\t")
    for index, row in df.iterrows():
        # ignore non-snp/indel variants
        if not pd.isnull(row['representative_transcript2']):
            continue
        if "-" in row['variant']:
            continue
        if 'EXPRESSION' in row['variant']:
            continue
        if 'fusion' in row['variant'].lower():
            continue
        if 'AMPLIFICATION' in row['variant']:
            continue
        if 'PHOSPHORYLATION' in row['variant']:
            continue
        if 'REARRANGEMENT' in row['variant']:
            continue
        if 'METHYLATION' in row['variant']:
            continue
        if 'TRANSLOCATION' in row['variant']:
            continue
        if 'SERUM' in row['variant']:
            continue
        if 'HOMOZYGOSITY' in row['variant']:
            continue
        if 'Alu insertion' in row['variant']:
            continue
        if 'MISLOCALIZATION' in row['variant']:
            continue
        if 'ALTER' in row['variant']:
            continue
        if 'VARIATION' in row['variant']:
            continue
        if 'WILD' in row['variant']:
            continue

        # match via genomic coordinate
        if row['drugs'] is np.nan:
                drug = ""
        else:
            drug = row['drugs']

        evidence_statement = str(row['evidence_statement'].encode("utf-8"))
        if not pd.isnull(row['reference_bases']):
            ikey = (row['chromosome'], int(row['start']), row['reference_bases'], row['variant_bases'])
            civic_variants[ikey] = (row['variant'], row['disease'], drug, row['evidence_level'], 
                                    evidence_statement, row['variant_origin'], 
                                    row['citation_id'], row['citation'])
            continue

        # rest of variant matches via HGVS momenclature
        m1 = re.search(r'[A-Z]\d+[A-Z]', row['variant'])
        m2 = re.search(r'[A-Z]\d+', row['variant'])
        m3 = re.search(r'([cC]\.(\d+)[\+-]\d+[A-Z]>[A-Z])', row['variant'])
        m4 = re.search(r'EXON (\d+) MUTATION', row['variant'])
        m5 = re.search(r'EXON (\d+)-(\d+) MUTATION', row['variant'])
        if m1 or m2:
            civic_phgvs[row['gene']][row['variant']] = (row['variant'], row['disease'], drug, row['evidence_level'], 
                                                        evidence_statement, row['variant_origin'], 
                                                        row['citation_id'], row['citation'])
        elif m3:
            variant = "c." + m3.group(0)[2:]
            civic_chgvs[row['gene']][variant] = (row['variant'], row['disease'], drug, row['evidence_level'], 
                                                 evidence_statement, row['variant_origin'], 
                                                 row['citation_id'], row['citation'])
        elif m4:
            civic_exon[row['gene']][m4.group(1)] = (row['variant'], row['disease'], drug, row['evidence_level'], 
                                                    evidence_statement, row['variant_origin'], 
                                                    row['citation_id'], row['citation'])
        elif m5:
            civic_exon[row['gene']][m5.group(1)] = (row['variant'], row['disease'], drug, row['evidence_level'], 
                                                    evidence_statement, row['variant_origin'], 
                                                    row['citation_id'], row['citation'])
            civic_exon[row['gene']][m5.group(2)] = (row['variant'], row['disease'], drug, row['evidence_level'], 
                                                    evidence_statement, row['variant_origin'], 
                                                    row['citation_id'], row['citation'])

    return civic_variants, civic_phgvs, civic_chgvs, civic_exon


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="VCF file that has HGVS names annotated")
    parser.add_argument("civic", help="CIVic database (ClinicalEvidenceSummaries)")
    parser.add_argument("-o", "--output", help="output filename (optional)")
    args = parser.parse_args()

    ifname = os.path.abspath(args.input)
    if not args.output:
        ofname = os.path.join(os.path.dirname(ifname), 
                              os.path.basename(ifname).split(".")[0] + ".civic_anno.vcf")
    else:
        ofname = args.output
    
    main(ifname, ofname, args.civic)



