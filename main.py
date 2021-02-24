import re
import matplotlib.pyplot as plt
from typing import Tuple
import numpy as np


def get_sequences(file_path: str) -> dict:
    genes = {}
    current_seq_list = []

    with open(file_path, 'r') as f:
        current_header = next(f)
        for raw_line in f:
            # Removing trailing newlines.
            line = raw_line.rstrip()
            # Checking if the current line is a header.
            if line[0] == '>':
                # Turning sequence list into one string.
                current_seq = ''.join(current_seq_list)
                try:
                    search = re.search(r'gene=[A-Za-z0-9]*',
                                       current_header)
                    gene_name = search.group().lower()
                except AttributeError:
                    gene_name = current_header

                genes[gene_name] = current_seq

                # Setting the current line as the new current
                # header.
                current_header = line
                # Resetting the current sequence list because
                # the previous sequence has ended.
                current_seq_list = []
            else:
                # Adding the current line to the sequence list.
                current_seq_list.append(line)

        # Turning sequence list into one string.
        current_seq = ''.join(current_seq_list)
        try:
            search = re.search(r'gene=[A-Za-z0-9]*', current_header)
            gene_name = search.group().lower()
        except AttributeError:
            gene_name = current_header

        genes[gene_name] = current_seq

        return genes


def separate(genes: dict) -> Tuple[dict, dict]:
    surface_protein = {}
    other = genes.copy()
    other.pop('gene=nef')
    surface_protein['gene=nef'] = genes['gene=nef']
    return surface_protein, other


order = ['A ', 'R', 'N', 'D', 'C', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
         'F', 'P', 'S', 'T', 'Trp', 'Y', 'V', 'Start', 'Stop']

aa3 = {"A ": ["GCT", "GCC", "GCA", "GCG"],
       "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
       "N": ["AAT", "AAC"], "D": ["GAT", "GAC"],
       "C": ["TGT", "TGC"], "E": ["CAA", "CAG", "GAA", "GAG"],
       "G": ["GGT", "GGC", "GGA", "GGG"],
       "H": ["CAT", "CAC"], "I": ["ATT", "ATC", "ATA"],
       "L": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
       "K": ["AAA", "AAG"], "M": ["ATG"], "F": ["TTT", "TTC"],
       "P": ["CCT", "CCC", "CCA", "CCG"],
       "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
       "T": ["ACT", "ACC", "ACA", "ACG"], "Trp": ["TGG"],
       "Y": ["TAT", "TAC"], "V": ["GTT", "GTC", "GTA", "GTG"],
       "Start": ["ATG", "CTG", "TTG", "GTG", "ATT"],
       "Stop": ["TAG", "TGA", "TAA"]}


def codon_counter(seq: str) -> dict:
    codon_count = {}
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        try:
            codon_count[codon] += 1
        except KeyError:
            codon_count[codon] = 1
    return codon_count


def sort_by_aa(codon_count: dict) -> dict:
    aa_count = {}
    for i in order:
        aa_count[i] = []
        for codon in aa3[i]:
            try:
                aa_count[i].append(codon_count[codon])
            except KeyError:
                aa_count[i].append(0)
    return aa_count


def plot(aa_count: dict) -> None:
    x = np.arange(64, step=3.2)
    for counter, amino_acid in zip(x, order):
        codon_counts = aa_count[amino_acid]
        name_list = aa3[amino_acid]
        pos = np.linspace(counter, counter+len(codon_counts)*0.4,
                          len(codon_counts))
        plt.bar(pos, codon_counts, 0.4)

    plt.xticks(x, order)
    plt.gcf().set_size_inches([12.8, 4.8])
    plt.show()


if __name__ == '__main__':
    genes = get_sequences('fasta/virus/siv seq mRNA.fasta')
    surface_protein, other = separate(genes)
    plot(sort_by_aa(codon_counter(surface_protein['gene=nef'])))


