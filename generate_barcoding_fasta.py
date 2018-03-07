from DNASkittleUtils.Contigs import read_contigs, write_contigs_to_file, Contig
from os.path import splitext
import sys

def informative_bases(seq, informative_positions):
    return ''.join([seq[i] for i in informative_positions])


def barcoding_sequences(input_fasta_path):
    consensus_contigs = read_contigs(input_fasta_path)
    con = consensus_contigs[0].seq
    informative_positions = [i for i, c in enumerate(con) if c not in 'ACGT']
    reduced_contigs = []
    for species in consensus_contigs[1:]:  # skip consensus
        seq = informative_bases(species.seq, informative_positions)
        reduced_contigs.append(Contig(species.name + '__reduced', seq))
    out_path = splitext(input_fasta_path)[0] + '__reduced_consensus.fa'
    write_contigs_to_file(out_path, reduced_contigs)
    print("Output reduced fasta file:", out_path)

    table_name = splitext(out_path)[0] + '__positions.csv'
    with open(table_name, 'w') as outfile:
        outfile.write('List of ambiguous bases output (0-index based)\n')
        for x in informative_positions:
            outfile.write("%i,%s\n" % (x, con[x]))
    print("Output ambiguous position table:", table_name)


if __name__ == '__main__':
    # consensus_file = r"D:\josiah\Documents\Research\Colleagues\Will Plumb\Barcoding-Fraxinus\data2\consensus align.fa"
    if len(sys.argv) < 2:
        print("Usage: ", sys.argv[0], '<genes_aligned_to_consensus.fa>')
        print("Your fasta file should contain the first entry as Consensus sequence with ambiguous bases")
        print("at the positions where the species sequences differ.  Each subsequent entry should be aligned")
        print("to the conensus in multiple sequence alignment.")
    else:
        barcoding_sequences(sys.argv[1])