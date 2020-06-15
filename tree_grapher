# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:03:01 2020

@author: Maia
"""

import pydot
from PIL import Image
from save_T import T

G = pydot.Dot(graph_type='digraph')

ROOT = list(T.keys())[0]

def tree_crawler(t,i,G):
    """Plots scheme of a tree.

    Keyword arguments:
    t -- Tree
    i -- Parameter identifier (integer)
    G -- graphic object
    """
    for e in list(t[i].keys())[2:]:
        node = pydot.Node(e,shape='circle')
        G.add_node(node)
        edge = pydot.Edge(i,e)
        G.add_edge(edge)
        if len(t[i][e]) > 0:
            tree_crawler(t[i],e,G)  # descenso recursivo a los infiernos...
    return G

G = tree_crawler(T,ROOT,G)


G.write_png('Tree.png')
im = Image.open('Tree.png')
im.show()
