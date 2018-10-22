from utils import *
from ete3 import NCBITaxa, TextFace, TreeStyle, faces
import argparse
import os
from collections import Counter

parser = argparse.ArgumentParser(description="Plot a tree with sequence content information "
                                             "from the Refseq release for a given list of taxonomy identifiers."
                                 )
parser.add_argument('-t',
                    dest='taxid',
                    help="A unique taxid for which all descendant information will be plotted. "
                         "e.g 10239 for viruses, 9604 for Hominidae (great apes)",
                    required=True
                    )
parser.add_argument('-r',
                    dest='render',
                    help='Path to file the image is to be written. Valid extensions are .svg, .pdf, .png '
                         'and are detected automatically',
                    required=False
                    )
parser.add_argument('-db',
                    dest='taxonomy_dir',
                    help='Dir with the instantiated NCBI taxonomy DB required from ete3',
                    required=False
                    )
parser.add_argument('-j',
                    dest='json_fp',
                    help='JSON file with the mapping of taxonomy ids to accessions',
                    required=True
                    )
parser.add_argument('-prune',
                    dest='prune_level',
                    type=str,
                    help='Level on which the tree will be pruned. (e.g. family)'
                         'Make sure this is at least a level deeper than the given taxid',
                    required=False)
parser.add_argument('-keep_accessions',
                    dest="keep_accessions",
                    type=str,
                    help="Comma separated list of accession prefixes to calculate seq_content for "
                         "e.g 'NC_, NZ_'",
                    required=False)
parser.add_argument('-prune_zeros',
                    dest="prune_zeros",
                    help="If this option is supplied along with 'prune', "
                         "nodes with 0 seq_content are removed during pruning",
                    required=False,
                    action='store_true'
                    )
parser.add_argument('-write_to_file',
                    dest="write",
                    help="Optionally write the tree in newick format",
                    required=False,
                    action='store_true'
                    )
parser.add_argument('-plot_pies',
                    dest='plot_pies',
                    required=False,
                    help="(EXPERIMENTAL) Plot piecharts below each node, representing the percentages of "
                         "the counts for each accession type that contributes to the total",
                    action='store_true')


def get_seq_content_from_dict(node_id, seq_map, keep_accessions=None):
    """
    Given a node id and a dictionary containing sequence information
    calculate the total amount of sequence content for the node.
    If the node is not in the dictionary, return 0
    :param node_id: A node id to get the sequence content
    :param seq_map: A dictionary of the form {node_id: {accession : seq_size, ... }, ...}
    :param keep_accessions: A tuple of accession prefixes, used for filtering. e.g ('NC_','NZ_CM').
                            A single accession should be passed as a string
    :return: The sum of sequence content or 0
    """
    seq_content = 0
    int_id = int(node_id)
    if int(node_id) in seq_map.keys():
        seqs = seq_map[int_id]
        # If no accessions to filter for is provided, sum the seqs
        if keep_accessions is None:
            seq_content += sum(seqs.values())
        else:
            for seq in seqs:
                if seq.startswith(keep_accessions):  # basic check?
                    seq_content += seqs[seq]
    return seq_content


def sum_node_sequences(n, seq_map, *keep_accessions):
    """
    Recursively calculate the seq_content of each node
    :param n: A node of the tree
    :param seq_map: A dictionary of the form {node_id: {accession : seq_size, ... }, ...}
    :param keep_accessions: A tuple of accession prefixes, used for filtering. e.g ('NC_','NZ_CM').
                            A single accession should be passed as a string
    :return: The total sequence content of a node, summed over its children and from the seq_map
    """

    if n.is_leaf():  # Base condition
        n.seq_content = get_seq_content_from_dict(n.name, seq_map, *keep_accessions)
        return n.seq_content
    else:  # Recurse to children
        for child in n.get_children():
            n.seq_content += sum_node_sequences(child, seq_map, *keep_accessions)

        n.seq_content = sum(c.seq_content for c in n.get_children())
        # In case the taxid of the node is not a leaf but contains sequence information,
        # update the value
        n.seq_content += get_seq_content_from_dict(n.name, seq_map, *keep_accessions)

    return n.seq_content


def count_node_accession_types(n, seq_map):
    """
    Count the number of accession types for a given node
    :param n: Input node
    :param seq_map: A dictionary that contains the sequence information
    :return: A dictionary of the form {acc_type: count,...}
    """
    node_id = int(n.name)
    acc_types = []
    if node_id in seq_map:
        acc_list = list(seq_map[node_id].keys())
        for acc in acc_list:
            fields = acc.split('_')
            if fields[0] == 'NZ':
                if fields[1].startswith('CP'):
                    acc_types.append('NZ_CP')
                elif fields[1].startswith('CM'):
                    acc_types.append('NZ_CM')
                else:
                    acc_types.append('NZ')
            else:
                acc_types.append(fields[0])

        acc_types_counter = dict(Counter(acc_types))
        return acc_types_counter
    else:
        return {}


def update_parent_accession_counter(n):
    parent = n.up
    if parent is not None:
        for acc_type in n.acc_counter.keys():
            if acc_type in parent.acc_counter.keys():
                parent.acc_counter[acc_type] += n.acc_counter[acc_type]
            elif (acc_type not in parent.acc_counter.keys()) and (acc_type in n.acc_counter.keys()):
                parent.acc_counter[acc_type] = n.acc_counter[acc_type]


def add_piechart_to_node(node):
    color_map = {"NC": "cyan",
                 "AC": "aquamarine",
                 "NZ": "black",
                 "NZ_CP": "gray",
                 "NZ_CM": "lightgray",
                 "NW": "antiquewhite",
                 "NG": "gold",
                 "NT": "teal"}

    acc_types = list(node.acc_counter.keys())
    acc_types_counts = list(node.acc_counter.values())
    acc_sum = sum(acc_types_counts)
    acc_percentages = [(v / acc_sum) * 100 for v in acc_types_counts]
    acc_colors = [color_map[acc_type] for acc_type in acc_types]
    node.add_face(faces.PieChartFace(acc_percentages, width=100, height=100, colors=acc_colors),
                  column=0, position="branch-bottom")


def nodes_to_prune(tree, level, prune_zeros=False):
    """
    Get the nodes that will be kept with pruning. Nodes at the same level that have a different rank are also returned
    :param tree: An ete3.TreeNode() instance
    :param level: A string for the taxonomy level to be pruned. With the NCBI taxonomy
    :param prune_zeros: Prune out nodes that have seq_content zero
    this value is stored in the rank attribute of the node
    :return: A list of nodes to keep
    """

    # Initialize a list of nodes to keep
    final_nodes = []
    for node in tree.traverse():
        if node.rank == level:
            if prune_zeros is True:
                if node.seq_content > 0:
                    final_nodes += [node.name]
            else:
                final_nodes += [node.name]
            # For some nodes, their sister nodes might not have been annotated as the specified level/rank,
            # yet their least common ancestor is the level's parent node
            sisters = node.get_sisters()
            for sister in sisters:
                if sister.rank != level:
                    if prune_zeros is True:
                        if sister.seq_content > 0:
                            final_nodes += [sister.name]
                    else:
                        final_nodes += [sister.name]

    return final_nodes


def is_magnitude(integer):
    """
    Get the magnitude of the given integer
    :param integer:
    :return: "Giga" for billions, "Mega" for millions, "Kilo" for thousands
    """
    if integer / 10 ** 9 >= 1:
        return "Giga"
    elif 1000 > integer / 10 ** 6 >= 1:
        return "Mega"
    else:
        return "Kilo"


def layout(node):
    """
    Annotate the node based on the seq_content feature. Billons are suffixed with Gbp,
    millions are suffixed with Mbp and thousands with "Kbp"
    :param node: A node
    :return:
    """
    # Display scientific names
    faces.add_face_to_node(TextFace(node.sci_name), node, 0)

    # Display sequence information
    if is_magnitude(node.seq_content) == "Giga":
        label = str(round(node.seq_content / 10 ** 9, 2)) + "Gbp"
    elif is_magnitude(node.seq_content) == "Mega":
        label = str(round(node.seq_content / 10 ** 6, 2)) + "Mbp"
    else:
        label = str(round(node.seq_content / 10 ** 3, 2)) + "Kbp"

    faces.add_face_to_node(TextFace(label), node, column=0, position="branch-top")


def main():
    args = parser.parse_args()

    if args.taxonomy_dir:
        ncbi = NCBITaxa(dbfile=os.path.join(args.taxonomy_dir, 'taxa.sqlite'))
    else:
        ncbi = NCBITaxa()

    taxids = [int(taxid.strip()) for taxid in args.taxid.split(',')]

    print("Making list of descendants for taxon: {}".format(taxids))
    descendants = create_full_taxa_list(taxids, ncbi, include_parent=True)

    print("Loading sequence content information")
    refseq_mapping = get_acc2taxid_map(args.json_fp)

    print("Filtering to selected taxa")
    refseq_mapping2 = dict((l, refseq_mapping[str(l)]) for l in descendants
                           if str(l) in refseq_mapping)

    # Remove the original dict to free up some memory
    del refseq_mapping

    kept_descendants = list(refseq_mapping2.keys())

    # get the tree topology
    tree = ncbi.get_topology(kept_descendants)

    # TO DO
    # Multiple traversals of the tree - Maybe use get_cached_content method?

    print("Adding seq_content attribute to all nodes")
    for node in tree.traverse():
        node.add_features(seq_content=0)
        if args.plot_pies:
            node.add_features(acc_counter=count_node_accession_types(node, refseq_mapping2))

    print("Calculating seq_content for nodes")
    # Update the seq_content recursively
    for node in tree.traverse():
        if args.keep_accessions:
            # Reformat provided string to accessions
            prefixes = [prefix.strip() for prefix in args.keep_accessions.split(',')]
            if len(prefixes) == 1:
                # When only one accession, keep its string
                keep_accessions = prefixes[0]
            else:
                # Turn them into a tuple
                keep_accessions = tuple(prefixes)
            node.seq_content = sum_node_sequences(node, refseq_mapping2, keep_accessions)
        else:
        #     # If no string is provided, use None
        #     # TO DO: better use optional argument for the functions
            node.seq_content = sum_node_sequences(node, refseq_mapping2, args.keep_accessions)

    if args.plot_pies:
        for node in tree.traverse(strategy='postorder'):
            update_parent_accession_counter(node)
            add_piechart_to_node(node)

    if args.prune_level:
        print("Pruning tree")
        final_nodes = nodes_to_prune(tree, args.prune_level, prune_zeros=args.prune_zeros)
        tree.prune(final_nodes)

    print("Plotting")
    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = layout

    # # This highlights node 9604
    # H = tree&"9604"
    # H.set_style(NodeStyle())
    # H.img_style["bgcolor"] = "indianred"

    # Render or just show the image
    if args.render:
        tree.render(args.render, dpi=300, tree_style=ts)
    else:
        tree.show(tree_style=ts)

    # Write a newick tree
    if args.write:
        dirname, outbase = os.path.split(args.render)
        out_name = str(outbase.split('.')[0]) + '.nw'
        outfile = os.path.join(dirname, out_name)
        tree.write(features=["seq_content"], outfile=outfile, format=1, format_root_node=True)

    if args.plot_pies:
        for node in tree.traverse():
            print(node.sci_name, node.acc_counter)
    # tree.show(tree_style=ts)


if __name__ == '__main__':
    main()
