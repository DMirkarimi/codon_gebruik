import re


def get_sequences(file_path):
    nucleotide_dict = {}
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

                nucleotide_dict[gene_name] = current_seq

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

        nucleotide_dict[gene_name] = current_seq

        # Setting the current line as the new current
        # header.
        current_header = line
        # Resetting the current sequence list because
        # the previous sequence has ended.
        current_seq_list = []
        return nucleotide_dict


def seperate(nucleotide_dict):
    surface_protein = {}
    other = nucleotide_dict.copy()
    other.pop('gene=nef')
    surface_protein['gene=nef'] = nucleotide_dict['gene=nef']
    return surface_protein, other


order = ['Ala ', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'GlT', 'Gly', 'His',
         'Ile', 'LeT', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp',
         'Tyr', 'Val', 'Start', 'Stop']

aa3 = {"Ala ": ["GCT", "GCC", "GCA", "GCG"],
       "Arg": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
       "Asn": ["AAT", "AAC"], "Asp": ["GAT", "GAC"],
       "Cys": ["TGT", "TGC"], "Gln": ["CAA", "CAG"],
       "GlT": ["GAA", "GAG"], "Gly": ["GGT", "GGC", "GGA", "GGG"],
       "His": ["CAT", "CAC"], "Ile": ["ATT", "ATC", "ATA"],
       "LeT": ["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"],
       "Lys": ["AAA", "AAG"], "Met": ["ATG"], "Phe": ["TTT", "TTC"],
       "Pro": ["CCT", "CCC", "CCA", "CCG"],
       "Ser": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
       "Thr": ["ACT", "ACC", "ACA", "ACG"], "Trp": ["TGG"],
       "Tyr": ["TAT", "TAC"], "Val": ["GTT", "GTC", "GTA", "GTG"],
       "Start": ["ATG", "CTG", "TTG", "GTG", "ATT"],
       "Stop": ["TAG", "TGA", "TAA"]}


def codon_counter(seq):
    codon_count = {}
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        try:
            codon_count[codon] += 1
        except KeyError:
            codon_count[codon] = 1
    return codon_count


def sort_by_AA(codon_count):
    AA_count = {}
    for i in order:
        AA_count[i] = []
        for codon in aa3[i]:
            try:
                AA_count[i].append(codon_count[codon])
            except KeyError:
                AA_count[i].append(0)
    print(AA_count)



if __name__ == '__main__':

    seqs = get_sequences('fasta/virus/siv seq mRNA.fasta')
    surface_protein, other = seperate(seqs)
    sort_by_AA(codon_counter(surface_protein['gene=nef']))


