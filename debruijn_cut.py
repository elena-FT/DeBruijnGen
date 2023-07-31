# debruijn_cut.py
import sys

def load_graph_from_file(f):
    """
    Load the graph from the input file.

    Parameters:
        f (str): Filepath of the input file.

    Returns:
        dict: Dictionary containing the merged nodes information.
        int: Value of k (size of k-mers).
    """
    merged_nodes = {}
    with open(f, 'r') as file:
        num_nodes, k = map(int, file.readline().split())
        for _ in range(num_nodes):
            node_info = file.readline().split(' ', 4)
            node_info = [info.rstrip() for info in node_info]
            node_id = int(node_info[0])
            sequence = node_info[1]
            multiplicity = int(node_info[2])
            predecessors = [int(pred) for pred in node_info[3].split(',') if pred != ''] if len(node_info) > 3 else []   # Handle nodes with no predecessors
            successors = [int(succ) for succ in node_info[4].split(',') if succ != ''] if len(node_info) > 4 else []  # Handle nodes with no successors
            merged_nodes[node_id] = [node_id, sequence, multiplicity, predecessors, successors]

    return merged_nodes, k

def write_graph_to_file(graph, k, outf):
    """
    Write the modified graph to the output file.

    Parameters:
        merged_nodes (dict): Dictionary containing the merged nodes information.
        k (int): Size of k-mers.
        outf (str): Filepath of the output file.
    """
    with open(outf, 'w') as outfile:
        outfile.write(f"{len(graph)} {k}\n")
        for node_id, node in graph.items():
            outfile.write(f"{node[0]} {node[1]} {node[2]} ")
            outfile.write(','.join(str(pred) for pred in node[3]))
            outfile.write(' ')
            outfile.write(','.join(str(succ) for succ in node[4]))
            outfile.write('\n')
            

def reverse_complement(sequence):
    """
    Get the reverse complement of a DNA sequence.

    Parameters:
        sequence (str): Input DNA sequence.

    Returns:
        str: Reverse complement of the input sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in reversed(sequence))

def is_tip(node, k, merged_nodes):
    """
    Check if a node is a tip in the graph.

    Parameters:
        node (list): List representing the node with its information.
        k (int): Size of k-mers.
        merged_nodes (dict): Dictionary containing the merged nodes information.

    Returns:
        bool: True if the node is a tip, False otherwise.
    """
    if len(node[3]) > 1:  # Node has more than one predecessor, not a tip
        return False
    visited = set()
    current = node
    while True:
        if len(current[4]) != 1:  # Node has more than one successor, not a tip
            return False
        visited.add(current[0])
        id_next_node = current[4][0]
        next_node = merged_nodes[id_next_node]
        if next_node[0] in visited:  # Cycle found, not a tip
            return False
        if len(visited) >= 2 * k:  # Tip size is greater than or equal to 2k, not a tip
            return False
        current = next_node
        if current[2] < max(merged_nodes[neighbor][2] for neighbor in current[3]):  # Tip multiplicity check
            break
    return True

def debruijn_cut(f, outf):
    """
    Perform the De Bruijn graph cut to remove tips and write the modified graph to an output file.

    Parameters:
        f (str): Filepath of the input file.
        outf (str): Filepath of the output file.
    """
    # Step 1: Load the graph from the input file
    merged_nodes, k = load_graph_from_file(f)

    # Step 2: Remove tips
    while True:
        removed_tips = False
        nodes_to_remove = set()

        for id_node, node in merged_nodes.items():
            if is_tip(node, k, merged_nodes):
                nodes_to_remove.add(id_node)
                pred = node[3]
                succ = node[4]
                if len(pred) != 0:
                    predecessor = merged_nodes[node[3][0]]
                    predecessor[4].remove(id_node)
                    if len(predecessor[4]) == 0:
                        nodes_to_remove.add(predecessor[0])
                if len(succ) != 0:
                    successor = merged_nodes[node[4][0]]
                    successor[3].remove(id_node)
                    if len(successor[3]) == 0:
                        nodes_to_remove.add(successor[0])
                removed_tips = True

        if not removed_tips:
            break

        # Remove the tips
        for node in nodes_to_remove:
            del merged_nodes[node]

    # Write the graph to output file
    write_graph_to_file(merged_nodes, k, outf)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 debruijn_cut.py input_file output_file")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        debruijn_cut(input_file, output_file)
