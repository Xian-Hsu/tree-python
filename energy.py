import numpy as np
#import sys
#import os
#sys.path.append(os.path.abspath("/data/xianhsu/others/computational_astrophysics/project"))
from classes import *
from constants import *
from evolve import *
# solve the energy and angular momentum in certain step
def cal_energy(particles, tree, ntot, time):
    # calculate all energy in certain step
    K = kinetic_energy(particles, ntot)
    U_dep = gravitational_potential_energy(particles, tree, ntot)
    #U_indep = 0.
    U_indep = indep_potential_energy(particles, ntot)
    U = U_dep + U_indep
    E = K + U
    return E, K, U

#def potential_energy(particles, ntot):
    # the position dependent potential energy with given function
#    left

def kinetic_energy(particles, ntot):
    E_k = 0.
    for i in range(ntot):
        E_k += 0.5 * particles[i].mas*(particles[i].vel[0]**2 + particles[i].vel[1]**2 + particles[i].vel[2]**2)
    return E_k

def gravitational_potential_energy(particles, tree, ntot):
    # consider the potential between particles, approximated with tree
    E_u = 0.
    for i in range(ntot):
        particle = np.array([particles[i].pos, particles[i].vel])
        #print(particle)
        U = np.zeros(1)
        get_potential(particle, tree, i, U)
        U[0] = U[0]*particles[i].mas
        #print(U[0], flush=True)
        E_u += U[0]
    E_u = E_u*0.5
    #print(E_u, flush=True)
    return E_u

def get_potential(particle, cell, pid, U):
    # get potential for single particle and the tree
    if check_theta(particle, cell, pid):
        if pid not in cell.par:
            add_potential(particle, cell, pid, U)
    else:
        for subcell in cell.sub:
            get_potential(particle, subcell, pid, U)

def add_potential(particle, cell, pid, U):
    # the particle mass is not considered here, this would times in whole potential function
    small_r = 1e-2
    r = np.sqrt((particle[0][0]-cell.apos[0])*(particle[0][0]-cell.apos[0])+
    (particle[0][1]-cell.apos[1])*(particle[0][1]-cell.apos[1])+
    (particle[0][2]-cell.apos[2])*(particle[0][2]-cell.apos[2]))
    #if (G*cell.mas/(r+small_r))!=0:
    #    U = U - G*cell.mas/(r+small_r)
    U[0] = U[0] -G*cell.mas/(r+small_r)
    #print("cell mass "+str(cell.mas)+" r "+ str(r+small_r))
    #print(-G*cell.mas/(r+small_r), flush=True)
    #print(E_u, flush=True)

def cal_z_angular_momentum(particles, n_tot):
    # compute the z direction angular momentum
    L = 0.
    for i in range(n_tot):
        L += np.cross([particles[i].pos[0],particles[i].pos[1]],[particles[i].vel[0],particles[i].vel[1]]) * particles[i].mas
    return L

def indep_potential_energy(particles, ntot):
    # add time independent potential like sun)
    small_r = 1e-2
    E_u = 0.
    for i in range(ntot):
        pos = np.array(particles[i].pos)
#        print(pos,flush=True)
        r = np.sqrt(pos[0]**2+pos[1]**2+pos[2]**2)
        E_u = E_u - G*particles[i].mas*Msun/(r+small_r)
    return E_u


