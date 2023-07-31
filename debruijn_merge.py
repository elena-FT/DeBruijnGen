# debruijn_merge.py
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
    Write the graph to the output file.

    Parameters:
        graph (dict): Dictionary containing the merged nodes information.
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

def get_unique_neighbors(node, merged_nodes):
    """
    Return a tuple (predecessor, successor) with the index of nodes that are the unique neighbors in a direction with the current node.

    Parameters:
        node (list): List representing the node with its information.
        merged_nodes (dict): Dictionary containing the merged nodes information.

    Returns:
        tuple: Tuple containing the index of the unique predecessor and the unique successor.
    """
    pred = node[3]
    succ = node[4]
    finalPred = -1
    finalSucc = -1
    
    if len(pred) == 1:
        succPred = merged_nodes[pred[0]][4] # list of predecessor's successors
        predPred = merged_nodes[pred[0]][3] # list of predecessor's predecessors
        if (len(succPred) == 1 and succPred[0] == node[0]) or (len(predPred) == 1 and predPred[0] == node[0]):
            finalPred = pred[0]      
    if len(succ) == 1:
        predSucc = merged_nodes[succ[0]][3] # list of successor's predecessors
        succSucc =  merged_nodes[succ[0]][4] # list of successor's successors
        if (len(predSucc) == 1 and predSucc[0] == node[0]) or (len(succSucc) == 1 and succSucc[0] == node[0]):
            finalSucc = succ[0]
            
    return (finalPred, finalSucc)

def merge_nodes_orientation(seq1, seq2, k):
    """
    Get the orientation and the right merged sequence for the new merged node.

    Parameters:
        seq1 (str): Sequence of the first node to merge.
        seq2 (str): Sequence of the second node to merge.
        k (int): Size of k-mers.

    Returns:
        tuple: Tuple containing the orientation and the merged sequence. Returns (None, None) if the sequences cannot be merged.
    """
    ln1 = len(seq1)
    ln2 = len(seq2)
    if seq1[ln1-k+1:ln1-1] == seq2[:k-2]:
        return ("FF", seq1 + seq2[k-1:])
    if reverse_complement(seq1)[ln1-k+1:ln1-1] == reverse_complement(seq2)[:k-2]:
        return ("RR", seq2 + seq1[k-1:])
    if seq1[ln1-k+1:ln1-1] == reverse_complement(seq2)[:k-2]:
        return ("FR", seq1 + reverse_complement(seq2)[k-1:])
    if reverse_complement(seq1)[ln1-k+1:ln1-1] == seq2[:k-2]:
        return ("RF",  reverse_complement(seq1) + seq2[k-1:])
    return (None, None)

def get_merged_nodes(node1, node2, new_merged_seq, orientation, merged_nodes, k):
    """
    Create the new merged node, with the right sequence and neighbors.

    Parameters:
        node1 (list): List representing the first node to merge with its information.
        node2 (list): List representing the second node to merge with its information.
        new_merged_seq (str): Sequence of the new merged node.
        orientation (str): Orientation of the merge ('FF', 'RR', 'FR', 'RF').
        merged_nodes (dict): Dictionary containing the merged nodes information.
        k (int): Size of k-mers.

    Returns:
        list: List representing the new merged node with its information.
    """
    new_merged_node = [[]] * 5
    new_merged_node[0] = node1[0]
    new_merged_node[1] = new_merged_seq
    new_merged_node[2] = node1[2] + node2[2]
    
    if orientation == 'FF':
        new_merged_node[3] = node1[3]
        new_merged_node[4] = node2[4]
    elif orientation == 'RR':
        new_merged_node[3] = node2[3]
        new_merged_node[4] = node1[4]
    elif orientation == 'RF':
        new_merged_node[3] = node1[4]
        new_merged_node[4] = node2[4]
    elif orientation == 'FR':
        new_merged_node[3] = node1[3]
        new_merged_node[4] = node2[3]
    
    if node1[0] in new_merged_node[3]:
        new_merged_node[3].remove(node1[0])
    if node2[0] in new_merged_node[3]:
        new_merged_node[3].remove(node2[0])
    if node1[0] in new_merged_node[4]:
        new_merged_node[4].remove(node1[0])
    if node2[0] in new_merged_node[4]:
        new_merged_node[4].remove(node2[0])
    
    return new_merged_node

def update_neighbors(merged_nodes, nodeDel, newNode_id):
    """
    Update neighbors when a node is deleted after a merge.

    Parameters:
        merged_nodes (dict): Dictionary containing the merged nodes information.
        nodeDel (list): List representing the node to be deleted with its information.
        newNode_id (int): Index of the new merged node replacing the deleted node.
    """
    preds = nodeDel[3]
    succs = nodeDel[4]
    
    for pred in preds:
        if pred != newNode_id:
            # Check whether the node is in the list of successors or predecessors of the current node.
            if nodeDel[0] in merged_nodes[pred][4]:
                idx = merged_nodes[pred][4].index(nodeDel[0])
                merged_nodes[pred][4][idx] = newNode_id
            elif nodeDel[0] in merged_nodes[pred][3]:
                idx = merged_nodes[pred][3].index(nodeDel[0])
                merged_nodes[pred][3][idx] = newNode_id
    
    for succ in succs:
        if succ != newNode_id:
            # Check whether the node is in the list of successors or predecessors of the current node.
            if nodeDel[0] in merged_nodes[succ][3]:
                idx = merged_nodes[succ][3].index(nodeDel[0])
                merged_nodes[succ][3][idx] = newNode_id
            elif nodeDel[0] in merged_nodes[succ][4]:
                idx = merged_nodes[succ][4].index(nodeDel[0])
                merged_nodes[succ][4][idx] = newNode_id

def debruijn_merge(f, outf):
    """
    Merge nodes representing reverse complements in the De Bruijn graph.

    Parameters:
        f (str): Filepath of the input file containing the De Bruijn graph.
        outf (str): Filepath of the output file to save the merged graph.
    """
    # Step 1: Load the graph from the input file
    merged_nodes, k = load_graph_from_file(f)
    
    # Step 2: Merge nodes representing reverse complements
    while True:
        merged_something = False

        for node_id, node in list(merged_nodes.items()):
            if node_id in merged_nodes:
                pred, succ = get_unique_neighbors(merged_nodes[node_id], merged_nodes)

                noChange = 1
                if pred != -1:
                    orientation, new_merged_seq = merge_nodes_orientation(merged_nodes[pred][1], merged_nodes[node_id][1], k)
                    if orientation:
                        new_merged_node = get_merged_nodes(merged_nodes[pred], merged_nodes[node_id], new_merged_seq, orientation, merged_nodes, k)
                        update_neighbors(merged_nodes, merged_nodes[node_id], pred)
                        del merged_nodes[node_id]
                        merged_nodes[pred] = new_merged_node
                        merged_something = True
                        noChange = 0

                if succ != -1 and noChange:
                    orientation, new_merged_seq = merge_nodes_orientation(merged_nodes[node_id][1], merged_nodes[succ][1], k)
                    if orientation:
                        new_merged_node = get_merged_nodes(merged_nodes[node_id], merged_nodes[succ], new_merged_seq, orientation, merged_nodes, k)
                        update_neighbors(merged_nodes, merged_nodes[succ], node_id)
                        del merged_nodes[succ]
                        merged_nodes[node_id] = new_merged_node
                        merged_something = True

        if not merged_something:
            break

    # Step 3: Write the graph to output file
    write_graph_to_file(merged_nodes, k, outf)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 debruijn_merge.py input_file output_file")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        debruijn_merge(input_file, output_file)
