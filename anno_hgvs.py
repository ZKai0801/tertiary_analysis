__doc__ = """
utilised LRG database to determine which transcript has the most clinic usage (or the longest transcript),
and subsequently determine the correct HGVS name for each variant
VEP was used in advance to annotate all possibilities.

This script will generate five novel infos in the vcf file: gene, transcript, exon, cHGVS and pHGVS.
It will also re-write vep annotated info (i.e. ANN) to retain only one element.

"""
__date__ = "2020/01/07"
__author__ = "Kai"
__version__ = "v1.2"


import pysam
import sys
import pandas as pd
import os
import argparse


def annotate_hgvs(infname, ofname, refflat):
    """
    Annotate HGVS name to each variant
    """
    clinic_transcript = read_clinic_transcript_database(args.lrg)
    with pysam.VariantFile(infname, "r") as vcfin:
        if "gene" not in vcfin.header.info:
            vcfin.header.info.add("gene", "1", "String", "Gene name")
        if "transcript" not in vcfin.header.info:
            vcfin.header.info.add("transcript", "1", "String", "Transcript with the most clinical impacts")
        if "exon" not in vcfin.header.info:
            vcfin.header.info.add("exon", "1", "String", "Variant occurred on which exon")
        if "cHGVS" not in vcfin.header.info:
            vcfin.header.info.add("cHGVS", "1", "String", "cHGVS nomenclature")
        if "pHGVS" not in vcfin.header.info:
            vcfin.header.info.add("pHGVS", "1", "String", "pHGVS nomenclature")
        with pysam.VariantFile(ofname, "w", header=vcfin.header) as vcfout:
            for record in vcfin:
                gene, transcript, exon, chgvs, phgvs = parse_vep_anno(record.info['ANN'], clinic_transcript, refflat)
                if chgvs:
                    chgvs = chgvs.split(":")[1]
                if phgvs:
                    phgvs = phgvs.split(":")[1].replace("%3D", "=")
                record.info['gene'] = gene
                record.info['transcript'] = transcript
                record.info['exon'] = exon
                record.info['cHGVS'] = chgvs
                record.info['pHGVS'] = phgvs
                record.info['ANN'] = tuple([i for i in record.info['ANN'] if transcript in i])
                vcfout.write(record)


def parse_vep_anno(vep_records, clinic_transcript, refflat):
    """
    Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|
    cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|
    DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|REFSEQ_MATCH|GENE_PHENO|HGVS_OFFSET|HGVSg
    """
    gene_symbols = list(set([i.split("|")[3] for i in vep_records]))

    # get the transcript with the most clinical impacts
    for each_gene in gene_symbols:
        if each_gene in clinic_transcript:
            transcript = clinic_transcript[each_gene][0].split(".")[0]
            for each_vep in vep_records:
                if transcript in each_vep:
                    if each_vep.split("|")[8]:
                        cds = "exon" + each_vep.split("|")[8].split("/")[0]
                    else:
                        cds = "intron" + each_vep.split("|")[9].split("/")[0]
                    return (each_gene, each_vep.split("|")[6], cds, 
                            each_vep.split("|")[10], each_vep.split("|")[11])
        else:
            continue
    
    if len(vep_records) == 1:
        if vep_records[0].split("|")[8]:
            cds = "exon" + vep_records[0].split("|")[8].split("/")[0]
        else:
            cds = "intron" + vep_records[0].split("|")[9].split("/")[0]
        return (gene_symbols[0], vep_records[0].split("|")[6], cds, 
                vep_records[0].split("|")[10], vep_records[0].split("|")[11])

    # get the longest transcript
    transcripts = [i.split("|")[6].split(".")[0] for i in vep_records]
    longest_transcript = determine_longest_transcript(transcripts, refflat)
    vep_record = [i for i in vep_records if longest_transcript in i][0]
    if vep_record.split("|")[8]:
        cds = "exon" + vep_record.split("|")[8].split("/")[0]
    else:
        cds = "intron" + vep_record.split("|")[9].split("/")[0]
    return (vep_record.split("|")[3], vep_record.split("|")[6], cds, 
            vep_record.split("|")[10], vep_record.split("|")[11])


def determine_longest_transcript(transcripts, refflat):
    """
    using refflat to determine which transcript should be used
    """
    longest_length = 0
    longest_transcript = ""
    with open(refflat, "r") as fh:
        for line in fh:
            line = line.strip().split()
            if line[1] in transcripts:
                length = int(line[5]) - int(line[4])
                if length > longest_length:
                    longest_length = length
                    longest_transcript = line[1]
    return longest_transcript


def read_clinic_transcript_database(lrg):
    """
    read the clinical_transcript.tsv to get the transcript with most
    clinical impacts or longest length for each gene 
    """
    clinic_transcript = {}
    df = pd.read_csv(lrg, sep="\t", header=None)
    for index, row in df.iterrows():
        clinic_transcript[row[0]] = (row[1], row[2])
    return clinic_transcript
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="VEP annoatated VCF file")
    parser.add_argument("lrg", help="clinical_transcripts.tsv")
    parser.add_argument("refflat", help="the refFlat.tsv")
    parser.add_argument("-o", "--ofname")
    args = parser.parse_args()

    infname = os.path.abspath(args.input)
    if not args.ofname:
        ofname = os.path.join(os.path.dirname(infname), os.path.basename(infname).split(".")[0] + ".step6.hgvs_anno.vcf")
    else:
        ofname = os.path.abspath(args.ofname)
    annotate_hgvs(infname, ofname, args.refflat)
