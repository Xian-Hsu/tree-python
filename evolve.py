import numpy as np
import matplotlib.pyplot as plt
#import os
#import sys
#sys.path.append(os.path.abspath("/data/xianhsu/others/computational_astrophysics/project"))
from classes import *
from tree_structure import *
from search_tree import *
from visualize import *
from constants import *

def rk4(particles, tree, pid, time, dt):
    # evolve with rk4 method, assume time indep. potential and without rebuild the tree
    # create a temporary particle with form [[pos], [vel], [acc]]
    par = np.array([particles[pid].pos, particles[pid].vel])
    k1 = get_diff(par, tree, time, pid)
    k2 = get_diff(par+0.5*dt*k1, tree, time+0.5*dt, pid)
    k3 = get_diff(par+0.5*dt*k2, tree, time+0.5*dt, pid)
    k4 = get_diff(par+dt*k3    , tree, time+dt    , pid)
    k  = 1/6*(k1+2*k2+2*k3+k4)
    for i in range(3):
        particles[pid].pos[i] += k[0][i]*dt
        particles[pid].vel[i] += k[1][i]*dt
    time += dt
    # add the value by the differential term
def mpi_rk4(particle, tree, pid, time, dt):
    # evolve with rk4 method, assume time indep. potential and without rebuild the tree
    # create a temporary particle with form [[pos], [vel], [acc]]
    par = np.array([particle.pos, particle.vel])
    k1 = get_diff(par, tree, time, pid)
    k2 = get_diff(par+0.5*dt*k1, tree, time+0.5*dt, pid)
    k3 = get_diff(par+0.5*dt*k2, tree, time+0.5*dt, pid)
    k4 = get_diff(par+dt*k3    , tree, time+dt    , pid)
    k  = 1/6*(k1+2*k2+2*k3+k4)
    for i in range(3):
        particle.pos[i] += k[0][i]*dt
        particle.vel[i] += k[1][i]*dt
    time += dt
    return particle
    # add the value by the differential term

def get_diff(particle, tree, time, pid):
    # update the velocity then consider all terms of the force sources
    # TODO double check this function
    diff = np.zeros([2,3])
    diff[0][:] = particle[1][:]
    get_gravity(particle, tree, pid, diff)
    get_potential(particle, diff)
    return diff

def get_gravity(particle, cell, pid, diff):
    # add the acceleration to particles, acc should be reset to 0 before it
    if check_theta(particle, cell, pid):
        if pid not in cell.par:
            add_gravity(particle, cell, pid, diff)
    else:
        for subcell in cell.sub:
            get_gravity(particle, subcell, pid, diff)

def check_theta(particle, cell, pid):
    # len / distance < theta is needed for calculate the force of the cell, if not, find subcells
    theta = 0.5
    check = False
    if (pid not in cell.par) and cell.len/np.sqrt((particle[0][0]-cell.apos[0])*(particle[0][0]-cell.apos[0])+
        (particle[0][1]-cell.apos[1])*(particle[0][1]-cell.apos[1])+
        (particle[0][2]-cell.apos[2])*(particle[0][2]-cell.apos[2])) < theta:
        check = True
    return check

def add_gravity(particle, cell, pid, diff):
    # directly add gravity from specific cell to the particle
    small_r = 1e4 # prevent infinite gravity
    rdif = np.array([particle[0][0]-cell.apos[0], particle[0][1]-cell.apos[1], particle[0][2]-cell.apos[2]])
    r = np.sqrt((particle[0][0]-cell.apos[0])*(particle[0][0]-cell.apos[0])+
    (particle[0][1]-cell.apos[1])*(particle[0][1]-cell.apos[1])+
    (particle[0][2]-cell.apos[2])*(particle[0][2]-cell.apos[2]))
    #print(rdif)
    for j in range(3):
        diff[1][j] =  diff[1][j] - G*cell.mas/(r+small_r)**2 * rdif[j]/(r+small_r)
#    if pid==1:
#        print("par, rdif, r, diff")
#        print(particle)
#        print(rdif)
#        print(r)
#        print(diff, flush=True)

def get_potential(particle, diff):
    small_r = 1e4
    r = np.sqrt(particle[0][0]**2+particle[0][1]**2+particle[0][2]**2)
    for j in range(3):
        diff[1][j] = diff[1][j] - G*Msun/(r+small_r)**2 * particle[0][j]/(r+small_r)

