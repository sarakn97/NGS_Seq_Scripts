#!/usr/bin/env python3
"""
Name: Sara K Nicholson
Title: Labeling HIV Sequence as X4 or R5 (X4 or R5-derived), returns dictionary of annotations
Date: February 24, 2025
Description: This program can be run on the command line by providing a fasta file with your sequences
and a length at which you would like to filter sequences at.
"""

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import subprocess

# First align sequences using subprocess module - subprocesses calls external packages
# (clustalw2 needs to be referenced or in path to be called here) - can use other MSAs
def clustal_align(infile, outfile):
    subprocess.run(["/home/sara/Downloads/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2", infile, outfile])
    # Convert Alignment file to Fasta

def convert_alignment(outfile, outfileFA)
    outFA = AlignIO.convert(outfile, "clustal", outfileFA, "fasta")
    return outFA

#  subprocess.run(["/home/sara/Downloads/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2", "-infile=/home/sara/GITHUB/Biopython_helper_scripts/sample_fasta.fa", "-outfile=/home/sara/GITHUB/Biopython_helper_scripts/aligned.aln"])
#    # Convert Alignment file to Fasta
#    AlignIO.convert("/home/sara/GITHUB/Biopython_helper_scripts/aligned.aln", "clustal", "/home/sara/GITHUB/Biopython_helper_scripts/aligned.fa", "fasta")


# View Aligned Sequences
alignment = AlignIO.read("/home/sara/GITHUB/Biopython_helper_scripts/aligned.fa", "fasta")
print(alignment)
for record in alignment:
    print(f"Sequence ID: {record.id}, Sequence: {record.seq}")

# second find divergence against a reference and have many if/then statements to conclude with label
# LOAD REFERENCES
refX4 = SeqRecord(
                Seq("TRPNNNTRKSIRIQRGPGRAFVTIGKI-GNMRQAHCNISRAKWNATLKQIASKLREQFGNNKTIIFKQSSGGDPEI"),
                id="NL43",
                name="REF_X4",
                description="sequence",
            )

refR5 = SeqRecord(
                Seq("TRPNNNTRKSINI--GPGRAFYTTGEIIGDIRQAHCNLSRAKWNDTLNKIVIKLREQFG-NKTIVFKHSSGGDPEI"),
                id="BAL",
                name="REF_R5",
                description="sequence",
            )

# prepare references so that each nucleotide is its own object in list
x4 = list(refX4.seq)
r5 = list(refR5.seq)

def list_sequences(recs, length):
    """ return list of sequences from SeqIO object of which are of certain length or greater """
    seq_rec = SeqIO.parse(recs, "fasta")
    fasta = []
    for record in seq_rec:
        if len(record.seq) >= length:  # remove sequences with <= 60 nucleotides
            record_final = SeqRecord(
                record.seq,
                id=record.id,
                name="ASV",
                description="sequence",
            )
            fasta.append(record_final)  # add sequence to list
    return fasta


def list_list_sequences(fa):
    """ return list of list of sequences """
    listed_seqs = []
    for record in fa:
        record = list(record.seq)
        listed_seqs.append(record)

    return listed_seqs


def annotate_seqs(seqs):
    """ score alignments to references and determine if X4 or R5 sequence based on score"""
    scoresx4 = []
    scoresr5 = []

    for seq in seqs:
        score = 0  # set base score zero
        for i in range(len(seq)):  # iterate through each position in sequence
            if seq[i] != x4[i]:  # score based on if it matches the x4 position (higher score = less similarity)
                score += 1
        scoresx4.append(score)

    for seq in seqs:
        score = 0
        for i in range(len(seq)):
            if seq[i] != r5[i]:
                score += 1
        scoresr5.append(score)

    annotation = []  # open list
    for i in range(len(scoresx4)):  # iterate
        if scoresx4[i] > scoresr5[i]:
            trop = "R5"  # if the x4 score is greater than the r5 score then that sequence is designated R5
        elif scoresx4[i] == scoresr5[i]:
            trop = "NA"  # if the scores are equal then the sequence tropism is indistinguishable
        else:
            trop = "X4"  # otherwise, the sequence is designated X4
        annotation.append(trop)  # add annotation designation to list

    print(scoresx4)
    print(scoresr5)
    print(annotation)
    return annotation


def get_seq_ids(fa):
    seqs =[]
    for record in fa:
        record = record.id
        seqs.append(record)

    return seqs


if __name__ == "__main__":
    # Initializations
    records = str(sys.argv[1])
    outfile = str(sys.argv[2])
    outfileFA = str(sys.argv[3])
    leng = int(sys.argv[4])

    clustal_align(records, outfile)
    aligned = convert_alignment(outfile, outfileFA)
    fasta = list_sequences(aligned, leng)
    listed_seqs = list_list_sequences(fasta)
    annotation = annotate_seqs(listed_seqs)
    seq_ids = get_seq_ids(fasta)

    # create dictionary using the list of sequence IDs and the designated annotations
    keys = seq_ids
    values = annotation

    # create dictionary using dict()
    seq_dict = dict(zip(keys, values))

    # print dictionary of annotations
    print(seq_dict)
