import os
import re
import sys

def sanitize_sequence(sequence):
    non_acgtn_chars = re.findall('[^ACGTN]', sequence)
    if non_acgtn_chars:
        print(f"Non-ACGTN characters found: {', '.join(set(non_acgtn_chars))}")
    return re.sub('[^ACGTN]', '', sequence)

def create_combined_fasta(input_dir, output_file):
    with open(output_file, 'w') as output_fasta:
        sequence_count = 0
        for subdir, _, files in os.walk(input_dir):
            for file in files:
                if file.endswith('.asm.bp.p_ctg.gfa.fa'):
                    file_path = os.path.join(subdir, file)
                    sequence_name = os.path.splitext(file)[0]
                    with open(file_path, 'r') as input_file:
                        sequence_data = input_file.read().strip()
                        sequences = sequence_data.split('>')
                        for sequence in sequences:
                            if sequence.strip():
                                sequence_count += 1
                                header, sequence_body = sequence.split('\n', 1)
                                output_fasta.write(f'>{sequence_name}_{sequence_count}\n')
                                output_fasta.write(f'{sequence_body}\n')

def read_fasta(file_path):
    """Read a FASTA file and return a dictionary of sequence IDs and sequences."""
    sequences = {}
    current_id = None
    with open(file_path, "r") as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                current_id = line[1:]
                sequences[current_id] = ""
            else:
                sequences[current_id] += line
    return sequences

def write_fasta(output_path, sequences):
    """Write a dictionary of sequences to a FASTA file."""
    with open(output_path, "w") as output_file:
        for seq_id, sequence in sequences.items():
            output_file.write(f">{seq_id}\n{sequence}\n")

# Get input directory, output file, HM7 sequences file, and output unused file paths from command line arguments
input_directory = sys.argv[1]
output_fasta_file = sys.argv[2]
hm7_sequences_file = sys.argv[3]
output_unused_file = sys.argv[4]

# Create combined fasta file
create_combined_fasta(input_directory, output_fasta_file)

# Read the sequences from HM7 sequences and the created combined fasta file
hm7_sequences = read_fasta(hm7_sequences_file)
combined_sequences = read_fasta(output_fasta_file)

# Find sequences that are in HM7 sequences but not in the combined fasta file
new_sequences = {}
for seq_id, sequence in hm7_sequences.items():
    if seq_id not in combined_sequences:
        new_sequences[seq_id] = sequence

# Write the new sequences to a new FASTA file
write_fasta(output_unused_file, new_sequences)
