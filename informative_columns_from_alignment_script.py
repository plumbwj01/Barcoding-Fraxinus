from sys import argv
from collections import Counter

# # Sparse Alignment columns

class Contig:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        
    def __repr__(self):
        return '< "%s" %i nucleotides>' % (self.name, len(self.seq))

def read_contigs(input_file_path):
    contigs = []
    current_name = ""
    seq_collection = []

    # Pre-read generates an array of contigs with labels and sequences
    with open(input_file_path, 'r') as streamFASTAFile:
        for read in streamFASTAFile.read().splitlines():
            if read == "":
                continue
            if read[0] == ">":
                # If we have sequence gathered and we run into a second (or more) block
                if len(seq_collection) > 0:
                    sequence = "".join(seq_collection)
                    seq_collection = []  # clear
                    contigs.append(Contig(current_name, sequence))
                current_name = read[1:]  # remove >
            else:
                # collects the sequence to be stored in the contig, constant time performance don't concat strings!
                seq_collection.append(read.upper())

    # add the last contig to the list
    sequence = "".join(seq_collection)
    contigs.append(Contig(current_name, sequence))
    return contigs


def collect_informative_column_information(input_file_path):
    species = read_contigs(input_file_path)
    informative_columns = {}
    consensus_sequence = []
    for col in range(len(species[0].seq)):
        letters = []
        for entry in species:
            letters.append(entry.seq[col])
        column_seq = ''.join(letters)
        consensusing = Counter(column_seq)
        consensus_sequence.append(consensusing.most_common()[0][0])
        if column_seq != letters[0] * len(species) and col > 100 and col < 1500:
            informative_columns[col] = column_seq
            print(column_seq, col + 1)
    species.append(Contig('Consensus', ''.join(consensus_sequence)))
    return informative_columns, species




# * Generate a fasta with informative columns
# * Majority vote consensus sequence, but it includes gaps
# * transpose?
# * CSV file write
def write_informative_columns_csv(input_file_path, output_name):
    informative_columns, species = collect_informative_column_information(input_file_path)
    with open(output_name, 'w') as csv_out:
        csv_out.write('Positions,' + ','.join([str(x + 1) for x in sorted(informative_columns.keys())]))
        csv_out.write('\n')
        for entry in species:
            csv_out.write(entry.name[:6] + ",")
            for col in range(len(species[0].seq)):
                if col in informative_columns:
                    csv_out.write(entry.seq[col] + ",")
            csv_out.write('\n')


if __name__ == '__main__':
    if len(argv) == 3:
        write_informative_columns_csv(argv[1], argv[2])
    else:
        print('Include input_path and output_path')
        print('Example: python informative_columns.py 9927_alignment.fasta 9927_informative_positions.csv')
        #
