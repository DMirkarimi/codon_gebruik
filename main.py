import re
import matplotlib.pyplot as plt
from typing import Tuple
import numpy as np
import pandas as pd
import ntpath


def get_sequences(file_path: str) -> dict:
    genes = {}
    current_seq_list = []

    with open(file_path, 'r') as f:
        current_header = next(f)
        for raw_line in f:
            # Removing trailing newlines.
            line = raw_line.rstrip()
            if line == '':
                continue
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


def separate(genes: dict) -> Tuple[list, list]:
    surface_protein = [genes['gene=nef']]
    genes.pop('gene=nef')
    other = list(genes.values())

    return surface_protein, other


order = ['A ', 'R', 'N', 'D', 'C', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
         'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'Start', 'Stop']

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
       "T": ["ACT", "ACC", "ACA", "ACG"], "W": ["TGG"],
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
    current_pos = 0
    pos_list = []
    # for count in range(0, 21):
    #     pos_list.append(current_pos)
    #     codon_counts = aa_count[order[count]]
    #
    #     pos = [current_pos+(0.5*i) for i in range(len(codon_counts))]
    #     plt.bar(pos, codon_counts, 0.4)
    #     current_pos += (len(codon_counts)+1) * 0.6
    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for count in range(0, 21):
        codon_counts = aa_count[order[count]]
        codon_names = aa3[order[count]]
        plt.bar(codon_names, codon_counts,
                label=order[count])

    plt.xticks(rotation=90)
    plt.legend(loc= 'lower center', ncol=10)
    # plt.xticks(pos_list, order)
    plt.gcf().set_size_inches([12, 5])
    name = ntpath.basename(file)
    plt.title(f'Codon Bias in {name}')
    plt.xlabel('Codons Per Aminozuur')
    plt.ylabel('Codon Bias in Procent')
    plt.tight_layout()
    plt.show()


def process_sequences(seqs: list) -> None:
    df = pd.DataFrame([codon_counter(seq) for seq in seqs])
    mean = dict(df.mean())
    for letter in order:
        aa = aa3[letter]
        total_count = 0
        for codon in aa:

            try:
                total_count += mean[codon]
            except:
                pass
        for codon in aa:
            try:
                x = mean[codon]/total_count *\
                                         100
                mean[codon] = x
            except:
                pass
    plot(sort_by_aa(mean))


if __name__ == '__main__':
    file = 'fasta/virus/hiv-1 seq mRNA.fasta'
    genes = get_sequences(file)
    surface_protein, other = separate(genes)
    process_sequences(surface_protein)
    process_sequences(other)

