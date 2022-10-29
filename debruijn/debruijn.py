#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import random
import statistics
import pickle
import networkx as nx
import matplotlib.pyplot as plt

random.seed(9001)

__author__ = "Mia Legras"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Mia Legras"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Mia Legras"
__email__ = "mialegras@gmail.com"
__status__ = "Developpement"

def isfile(path):
    '''
    Check if path is an existing file.

    Parameter
    ---------
    path : str
        Path to the file

    Return
    ------
    str
        Path if is an existing path to a file.
    '''
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    '''
    Retrieves the arguments of the program.

    Returns
    -------
    An object that contains the arguments.
    '''
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    '''
    Read sequence fasta from fastq file.

    Parameter
    ---------
    fastq_file : str
        Name of the fastq file.

    Return
    ------
    generator
        Generator of the sequence as string.
    '''
    isfile(fastq_file)
    with open(fastq_file, 'r') as fastq:
        for line in fastq:
            yield next(fastq).strip()
            next(fastq)
            next(fastq)


def cut_kmer(read, kmer_size):
    '''
    Cut a sequence into subsequence of the size kmer_size.

    Parameters
    ----------
    read : str
        Sequence to cut.
    kmer_size : int
        Size of kmer.

    Return
    ------
    generator
        Generator of list of subsequences.
    '''
    for i in range(len(read)):
        yield read[i:i+kmer_size]



def build_kmer_dict(fastq_file, kmer_size):
    '''
    Build a dictionary of kmer with their occurences.

    Parameters
    ----------
    fastq_file : str
        Name of the fastq file.
    kmer_size : int
        Size of kmer.

    Return
    ------
    dict
        A dictionary with kmer as key and their occurrence as value.
    '''
    dico = {}
    seq = "".join(list(read_fastq(fastq_file)))
    kmer = list(cut_kmer(seq, kmer_size))
    for each in kmer:
        if each in dico and len(each) == kmer_size:
            dico[each] += 1
        if each not in dico and len(each) == kmer_size:
            dico[each] = 1
    return dico


def build_graph(kmer_dict):
    '''
    Build a graph with a node for each kmer.

    Parameter
    ---------
    kmer_dict : dict
        Dictionary of kmer.

    Return
    ------
    graph
        A Digraph where nodes represent kmers.
    '''
    G = nx.DiGraph()
    for kmer in kmer_dict:
        G.add_edge(kmer[:-1], kmer[1:], weight=kmer_dict[kmer])
    return G


def remove_paths(graph, path_list,
                    delete_entry_node=True, delete_sink_node=True):
    '''
    Remove paths in graph.

    Parameters
    ----------
    graph : graph
        A Digraph where nodes represent kmers.
    path_list : list
        List of paths.
    delete_entry_node : boolean
        If True the entry node is removed.
    delete_sink_node : boolean
        If True the sink node is removed.

    Return
    ------
    graph
        A new Digraph where some paths were removed.
    '''
    for path in path_list:
        if delete_entry_node and delete_sink_node is not True:
            graph.remove_nodes_from(path[:-1])
        elif delete_entry_node is not True and delete_sink_node:
            graph.remove_nodes_from(path[1:])
        elif delete_entry_node is not True and delete_sink_node is not True:
            graph.remove_nodes_from(path[1:-1])
        else:
            graph.remove_nodes_from(path)
        return graph


def std(data):
    '''
    Compute standard-error on data.

    Parameter
    ---------
    data : list
        List of numerics values.

    Return
    ------
    list
        List of standard-error.'''
    return statistics.stdev(data)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    weight = std(weight_avg_list)
    length = std(path_length)
    if weight > 0:
        best_ind = weight_avg_list.index(max(weight_avg_list))
    elif weight == 0:
        if length > 0:
            best_ind = path_length.index(max(path_length))
        elif length == 0:
            best_ind = random.randint(0, len(path_list)-1)
    for i in range(len(path_list)):
        if path_list[i] != path_list[best_ind]:
            graph = remove_paths(graph, [path_list[i]],
                delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph, path):
    sub = graph.subgraph(path).edges(data=True)
    return statistics.mean([d["weight"] for (u, v, d) in sub])


def solve_bubble(graph, ancestor_node, descendant_node):
    if nx.has_path(graph, ancestor_node, descendant_node):
        paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
        list_weight = []
        list_len = []
        for path in paths:
            list_len.append(len(path))
            list_weight.append(path_average_weight(graph, path))
        if len(list_weight) > 1:
            return select_best_path(graph, paths, list_len, list_weight)
            #remove_paths(graph, path, ancestor_node, descendant_node)
    return graph


def simplify_bubbles(graph):
    bubble = False
    for node in graph.nodes:
        list_pred = list(graph.predecessors(node))
        if len(list_pred) > 1:
            for i, pred1 in enumerate(list_pred):
                for pred2 in list_pred[:1] + list_pred[i+1:]:
                    ancestor_node = nx.lowest_common_ancestor(graph,
                                                            pred1, pred2)
                    if ancestor_node is not None:
                        bubble = True
                        break
        if bubble:
            break
    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, ancestor_node, node))
    return graph


def solve_entry_tips(graph, starting_nodes):
    path_list = []
    pred_start = 0
    for node in graph:
        preds_list = list(graph.predecessors(node))
        if len(preds_list) > 1:
            for start in starting_nodes:
                path = list(nx.all_simple_paths(graph, start, node))[0]
                if path:
                    pred_start += 1
                    path_list.append(path)
                    if pred_start == 2 :
                        break
        if pred_start == 2:
            break
    # s'il y a plus de 2 starting nodes
    if pred_start == 2 :
        weight_list = []
        path_len = []
        for path in path_list:
            weight_list.append(path_average_weight(graph, path))
            path_len.append(len(path))

        graph = select_best_path(graph, path_list, path_len, weight_list,
                        True, False)
        starting_nodes = get_starting_nodes(graph)
        graph = solve_entry_tips(graph, starting_nodes)
    return graph


def solve_out_tips(graph, ending_nodes):
    pass


def get_starting_nodes(graph):
    nodes = list(graph.nodes())
    entry = []
    for node in nodes:
        if len(list(graph.predecessors(node))) == 0:
            entry.append(node)
    return entry


def get_sink_nodes(graph):
    nodes = list(graph.nodes())
    end = []
    for node in nodes:
        if len(list(graph.successors(node))) == 0:
            end.append(node)
    return end


def get_contigs(graph, starting_nodes, ending_nodes):
    contigs = []
    for start_nod in starting_nodes:
        for end_nod in ending_nodes:
            if nx.has_path(graph, start_nod, end_nod):
                for path in nx.all_simple_paths(graph, start_nod, end_nod):
                    seq = path[0]
                    for node in path[1:]:
                        seq+=node[-1]
                    contigs.append(tuple((seq, len(seq))))
    return contigs


def save_contigs(contigs_list, output_file):
    with open(output_file, 'w') as fasta:
        for i in range(len(contigs_list)):
            fasta.write(f'>contig_{i} len={contigs_list[i][1]}\n')
            fasta.write(fill(contigs_list[i][0]))
            fasta.write('\n')


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    seq = read_fastq(args.fastq_file)
    print(seq)

    #kmer = cut_kmer('TCAGAG', 3)
    
    dico = build_kmer_dict(args.fastq_file, 22)
    graph = build_graph(dico)
    #subax1 = plt.subplot(121)
    nx.draw(graph, with_labels=True, font_weight='bold')
    #subax2 = plt.subplot(122)
    nx.draw_shell(graph, with_labels=True, font_weight='bold')
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit 
    # graphe
    # Plot the graph
    # if args.graphimg_file:
    #     draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    # if args.graph_file:
    #     save_graph(graph, args.graph_file)


if __name__ == '__main__':
    main()
