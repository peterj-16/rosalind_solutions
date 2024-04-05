import re 

table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
        'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }

# read in fasta
def read_fasta(file_path):
    sequences = {}
    current_header = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_header = line[1:]
                sequences[current_header] = ''
            else:
                sequences[current_header] += line
    return sequences

file_path = "path_to_fasta"

# function to translate DNA to protein
def translate(seq):
    protein =''
    if len(seq)%3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein+= table[codon]
    return protein


# Find and translate genes in the original sequences
for header, seq in read_fasta(file_path).items():
    # Find all start codon positions in the original sequence
    start_codon_indices = [m.start() for m in re.finditer('ATG', seq)]

    for start_codon_index in start_codon_indices:
        # Find all occurrences of stop codons after the start codon in the original sequence
        stop_codon_pattern = re.compile('(TAA|TAG|TGA)')
        stop_codon_matches = stop_codon_pattern.finditer(seq, start_codon_index + 3)

        for stop_codon_match in stop_codon_matches:
            stop_codon_index = stop_codon_match.start()

            # Extract the gene sequence from the original sequence
            gene_sequence = seq[start_codon_index:stop_codon_index + 3]

            # Ensure the length of the gene sequence is a multiple of 3
            gene_sequence_length = len(gene_sequence)
            if gene_sequence_length % 3 == 0:
                # Translate the gene sequence into a protein sequence
                protein_sequence = translate(gene_sequence)

                # Print the results
                print(f"Header: {header}")
                print(f"Gene Sequence: {gene_sequence}")
                print(f"Protein Sequence: {protein_sequence}\n")

                # Stop processing this start codon after finding the first valid stop codon
                break

    # Find all start codon positions in the reverse complement
    reverse_complement = ''.join([{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}[base] for base in reversed(seq)])
    reverse_start_codon_indices = [m.start() for m in re.finditer('ATG', reverse_complement)]

    for reverse_start_codon_index in reverse_start_codon_indices:
        # Find all occurrences of stop codons after the start codon in the reverse complement sequence
        reverse_stop_codon_pattern = re.compile('(TAA|TAG|TGA)')
        reverse_stop_codon_matches = reverse_stop_codon_pattern.finditer(reverse_complement, reverse_start_codon_index + 3)

        for reverse_stop_codon_match in reverse_stop_codon_matches:
            reverse_stop_codon_index = reverse_stop_codon_match.start()

            # Extract the gene sequence from the reverse complement sequence
            gene_sequence = reverse_complement[reverse_start_codon_index:reverse_stop_codon_index + 3]

            # Ensure the length of the gene sequence is a multiple of 3
            gene_sequence_length = len(gene_sequence)
            if gene_sequence_length % 3 == 0:
                # Translate the gene sequence into a protein sequence
                protein_sequence = translate(gene_sequence)

                # Print the results
                print(f"Header: {header}")
                print(f"Gene Sequence: {gene_sequence}")
                print(f"Protein Sequence: {protein_sequence}\n")

                # Stop processing this start codon after finding the first valid stop codon
                break
