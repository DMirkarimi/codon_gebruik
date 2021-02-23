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

aa3 = {"Ala ": ["GCU" , "GCC" , "GCA" , "GCG"] ,
       "Arg": ["CGU" , "CGC" , "CGA" , "CGG" , "AGA" , "AGG"] ,
       "Asn": ["AAU" , "AAC"] ,
       "Asp": ["GAU" , "GAC"] ,
       "Cys": ["UGU" , "UGC"] ,
       "Gln": ["CAA" , "CAG"] ,
       "Glu": ["GAA" , "GAG"] ,
       "Gly": ["GGU" , "GGC" , "GGA" , "GGG"] ,
       "His": ["CAU" , "CAC"] ,
       "Ile": ["AUU" , "AUC" , "AUA"] ,
       "Leu": ["UUA" , "UUG" , "CUU" , "CUC" , "CUA" , "CUG"] ,
       "Lys": ["AAA" , "AAG"] ,
       "Met": ["AUG"] ,
       "Phe": ["UUU" , "UUC"] ,
       "Pro": ["CCU" , "CCC" , "CCA" , "CCG"] ,
       "Ser": ["UCU" , "UCC" , "UCA" , "UCG" , "AGU" ,"AGC"] ,
       "Thr": ["ACU" , "ACC" , "ACA" , "ACG"] ,
       "Trp": ["UGG"] ,
       "Tyr": ["UAU" , "UAC"] ,
       "Val": ["GUU" , "GUC" , "GUA" , "GUG"] ,
       "Start": ["AUG" , "CUG" , "UUG" , "GUG" , "AUU"] ,
       "Stop" : ["UAG" , "UGA" , "UAA"]
}


if __name__ == '__main__':
    x = get_sequences('fasta/virus/siv seq mRNA.fasta')
    print(seperate(x)[0])

