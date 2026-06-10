#!/usr/bin/env python3

from Bio import SeqIO
import os

print("Genome\tGenome_size_bp\tGenome_size_Mbp\tNum_contigs")

for fasta in open("genome_files.txt"):
    fasta = fasta.strip()

    genome = os.path.basename(fasta)
    genome = genome.replace(".fna", "")

    total_bp = 0
    contigs = 0

    for record in SeqIO.parse(fasta, "fasta"):
        total_bp += len(record.seq)
        contigs += 1

    print(
        f"{genome}\t{total_bp}\t{total_bp/1e6:.3f}\t{contigs}"
    )
