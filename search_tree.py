import numpy as np
#import os
#import sys
#sys.path.append(os.path.abspath("/data/xianhsu/others/computational_astrophysics/project"))
from classes import *

def tree_1dlist(tree):
    # search each end of the node the tree and get a 1d list with all cells
    cell_list = []
    cell_list = tree_list_append(tree, cell_list)
    return cell_list

def tree_list_append(cell, cell_list):
    # this appends the cells to the 1dlist
    if cell.sub != []:
        for subcell in cell.sub:
            cell_list = tree_list_append(subcell, cell_list)
    else:
        cell_list = np.append(cell_list, [cell])
    return cell_list

