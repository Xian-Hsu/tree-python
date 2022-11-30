import numpy as np
import random as rd
#import os
#import sys
#sys.path.append(os.path.abspath("/data/xianhsu/others/computational_astrophysics/project"))
from classes import *

# import parameter : crit_ratio, tot_mass, box_length, particle_list

#crit = 1.e-3 # (should come from parameter file) critical ratio for cells: n_crit = n_tot * crit
#n_crit = n_tot * crit
def tree_structure(particles, max_level, n_tot, box_length, n_crit):
    # get data points, then create tree structure for each step
    tree = Cell(np.zeros(3), box_length)
    tree.par = np.append(tree.par, range(n_tot))
    min_length = box_length / 2**(max_level)
    tree_node(tree, n_crit, particles, min_length)
    return tree

def tree_node(cell, n_crit, particles, min_length):
    #print(n_crit, cell.len, min_length)
    if (len(cell.par) > 1 and cell.len > min_length):
        # record the data of the cell then divide to 8 subcell, finally clear the particle list
        cell.mas = 0
        xpos = 0
        ypos = 0
        zpos = 0
        if len(cell.par) != 0:
            for ii in cell.par:
                i = int(ii)
                cell.mas += particles[i].mas
                xpos += particles[i].pos[0] * particles[i].mas
                ypos += particles[i].pos[1] * particles[i].mas
                zpos += particles[i].pos[2] * particles[i].mas
            cell.apos[0] = xpos / cell.mas
            cell.apos[1] = ypos / cell.mas
            cell.apos[2] = zpos / cell.mas
        # cut to 8 subcells with name cell.sub[0~7]
        # rule: - for 0 and + for 1 (left to right), ex. index 5 = sum(+-+)
        #print("before ", cell.pos, cell.len)
        cell.sub = np.append(cell.sub, [
        Cell(np.array([cell.pos[0]-cell.len/4, cell.pos[1]-cell.len/4, cell.pos[2]-cell.len/4]), cell.len/2),
        Cell(np.array([cell.pos[0]+cell.len/4, cell.pos[1]-cell.len/4, cell.pos[2]-cell.len/4]), cell.len/2),
        Cell(np.array([cell.pos[0]-cell.len/4, cell.pos[1]+cell.len/4, cell.pos[2]-cell.len/4]), cell.len/2),
        Cell(np.array([cell.pos[0]+cell.len/4, cell.pos[1]+cell.len/4, cell.pos[2]-cell.len/4]), cell.len/2),
        Cell(np.array([cell.pos[0]-cell.len/4, cell.pos[1]-cell.len/4, cell.pos[2]+cell.len/4]), cell.len/2),
        Cell(np.array([cell.pos[0]+cell.len/4, cell.pos[1]-cell.len/4, cell.pos[2]+cell.len/4]), cell.len/2),
        Cell(np.array([cell.pos[0]-cell.len/4, cell.pos[1]+cell.len/4, cell.pos[2]+cell.len/4]), cell.len/2),
        Cell(np.array([cell.pos[0]+cell.len/4, cell.pos[1]+cell.len/4, cell.pos[2]+cell.len/4]), cell.len/2)])
        #print("after ", cell.pos, cell.len)
        # fill particle index to new subcell
        for ii in cell.par:
            i = int(ii)
            if   particles[i].pos[0]<cell.pos[0] and particles[i].pos[1]<cell.pos[1] and particles[i].pos[2]<cell.pos[2]:
                cell.sub[0].par = np.append(cell.sub[0].par, [i])
            elif particles[i].pos[0]>=cell.pos[0] and particles[i].pos[1]<cell.pos[1] and particles[i].pos[2]<cell.pos[2]:
                cell.sub[1].par = np.append(cell.sub[1].par, [i])
            elif particles[i].pos[0]<cell.pos[0] and particles[i].pos[1]>=cell.pos[1] and particles[i].pos[2]<cell.pos[2]:
                cell.sub[2].par = np.append(cell.sub[2].par, [i])
            elif particles[i].pos[0]>=cell.pos[0] and particles[i].pos[1]>=cell.pos[1] and particles[i].pos[2]<cell.pos[2]:
                cell.sub[3].par = np.append(cell.sub[3].par, [i])
            elif particles[i].pos[0]<cell.pos[0] and particles[i].pos[1]<cell.pos[1] and particles[i].pos[2]>=cell.pos[2]:
                cell.sub[4].par = np.append(cell.sub[4].par, [i])
            elif particles[i].pos[0]>=cell.pos[0] and particles[i].pos[1]<cell.pos[1] and particles[i].pos[2]>=cell.pos[2]:
                cell.sub[5].par = np.append(cell.sub[5].par, [i])
            elif particles[i].pos[0]<cell.pos[0] and particles[i].pos[1]>=cell.pos[1] and particles[i].pos[2]>=cell.pos[2]:
                cell.sub[6].par = np.append(cell.sub[6].par, [i])
            elif particles[i].pos[0]>=cell.pos[0] and particles[i].pos[1]>=cell.pos[1] and particles[i].pos[2]>=cell.pos[2]:
                cell.sub[7].par = np.append(cell.sub[7].par, [i])
            else:
                print("error happens")
        cell.par = [] # clean original particle index
        for i in range(8):
            tree_node(cell.sub[i], n_crit, particles, min_length)
    else:
        # this is the end of one node. now sum over all mass and get the position
        cell.mas = 0
        xpos = 0
        ypos = 0
        zpos = 0
        if len(cell.par) != 0:
            for ii in cell.par:
                i = int(ii)
                cell.mas += particles[i].mas
                xpos += particles[i].pos[0] * particles[i].mas
                ypos += particles[i].pos[1] * particles[i].mas
                zpos += particles[i].pos[2] * particles[i].mas
            cell.apos[0] = xpos / cell.mas
            cell.apos[1] = ypos / cell.mas
            cell.apos[2] = zpos / cell.mas
        else:
            for i in range(3):
                cell.apos[i] = cell.pos[i]
        for i in range(3):
            cell.pos[i] -= cell.len/2
    return cell
# parts under are tests

#n_tot = int(1e4)
#n_crit = n_tot * crit
#box_length = 1000
#max_level = 8
#min_length = box_length / 2**(max_level)
#particles = []

#for i in range(n_tot):
#    particles = np.append(particles, [Particle(1, [rd.uniform(-500, 500), rd.uniform(-500,500), rd.uniform(-500,500)], np.zeros(3), np.zeros(3))])
#tree = tree_structure(particles, max_level)
#print(tree)

