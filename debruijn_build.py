# debruijn_build.py
import sys

def debruijn_build(f, k, outf):
    """
    Build the De Bruijn graph from a file containing DNA sequences.

    Parameters:
        f (str): Filepath of the input file containing DNA sequences.
        k (int): Size of k-mers.
        outf (str): Filepath of the output file to save the De Bruijn graph.
    """

    # Define a function to get the reverse complement of a DNA sequence.
    def reverse_complement(s):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join(complement[base] for base in reversed(s))

    # Initialize the nodes dictionary to store k-mers and their properties.
    nodes = {}

    # Read the input file containing DNA sequences.
    with open(f, 'r') as file:
        read_id = -1
        for line in file:
            line = line.strip()

            # Check if the line starts with '@' indicating the start of a new read.
            if line.startswith('@'):
                read_id += 1
            else:
                sequence = line

                # Iterate over the sequence to extract k-mers and build the graph.
                for i in range(len(sequence) - k + 1):
                    kmer = sequence[i:i + k]
                    complement = reverse_complement(kmer)

                    # If the reverse complement exists in nodes, use it as the kmer.
                    if complement in nodes:
                        kmer = complement

                    # If the kmer already exists in nodes, increase its multiplicity.
                    if kmer in nodes:
                        nodes[kmer][2] += 1
                    else:
                        # Create a new node for the kmer if it doesn't exist.
                        nodes[kmer] = [len(nodes), kmer, 1, [], []]

                    # Add predecessors and successors for the new node.
                    if i > 0:
                        predecessor = sequence[i - 1:i - 1 + k]
                        pred_complement = reverse_complement(predecessor)

                        # If the reverse complement exists in nodes, use it as the predecessor.
                        if pred_complement in nodes:
                            predecessor = pred_complement

                        # If the predecessor kmer exists in nodes and is not already a predecessor or successor of the current kmer, add it as a predecessor.
                        if predecessor in nodes and not nodes[predecessor][0] in nodes[kmer][3] and not nodes[predecessor][0] in nodes[kmer][4]:
                            if len(nodes[kmer][3]) != 0 and len(nodes[kmer][4]) == 0:
                                nodes[kmer][4].append(nodes[predecessor][0])
                            else:
                                nodes[kmer][3].append(nodes[predecessor][0])

                    if i < len(sequence) - k:
                        successor = sequence[i + 1:i + 1 + k]
                        succ_complement = reverse_complement(successor)

                        # If the reverse complement exists in nodes, use it as the successor.
                        if succ_complement in nodes:
                            successor = succ_complement

                        # If the successor kmer doesn't exist in nodes and its reverse complement also doesn't exist, add it as a successor.
                        elif not successor in nodes and not succ_complement in nodes:
                            if len(nodes[kmer][4]) != 0 and len(nodes[kmer][3]) == 0:
                                nodes[kmer][3].append(len(nodes))
                            else:
                                nodes[kmer][4].append(len(nodes))

                        # If the successor kmer exists in nodes and is not already a predecessor or successor of the current kmer, add it as a successor.
                        elif successor in nodes and not nodes[successor][0] in nodes[kmer][4] and not nodes[successor][0] in nodes[kmer][3]:
                            nodes[kmer][4].append(nodes[successor][0])

    # Write the graph to the output file.
    with open(outf, 'w') as outfile:
        outfile.write(f"{len(nodes)} {k}\n")
        for node_id, node in nodes.items():
            outfile.write(f"{node[0]} {node[1]} {node[2]} ")
            outfile.write(','.join(str(pred) for pred in node[3]))
            outfile.write(' ')
            outfile.write(','.join(str(succ) for succ in node[4]))
            outfile.write('\n')

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 debruijn_build.py input_file k output_file")
    else:
        input_file = sys.argv[1]
        k = int(sys.argv[2])
        output_file = sys.argv[3]
        debruijn_build(input_file, k, output_file)
