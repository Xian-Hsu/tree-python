import numpy as np
import math
import random as rd
from constants import *
#import os
#import sys
#sys.path.append(os.path.abspath("/data/xianhsu/others/computational_astrophysics/project"))
from classes import *
from mpi4py import MPI
comm = MPI.COMM_WORLD

def initial(n_tot, particles):
    # initialize the particles with Monte-Carlo like method, first try with pure random condition
    for i in range(n_tot):
        check = 1
        while(check):
            position = np.array([rd.uniform(-500, 500), rd.uniform(-500,500), rd.uniform(-500,500)])
            if (position[0]**2+position[1]**2+position[2]**2) <= 250000:
                check = 0
        velocity = np.zeros(3)
        #velocity = np.array([rd.uniform(-5e-5, 5e-5), rd.uniform(-5e-1,5e-1), rd.uniform(-5e-1,5e-1)])
        particles = np.append(particles, [Particle(1e10, position, velocity)])
    return particles

def initial_count():
    # create string array for plotting
    count = []
    for i in range(10):
        count = np.append(count, [str('000')+str(i)])
    for i in range(10,100):
        count = np.append(count, [str('00') +str(i)])
    for i in range(100, 1000):
        count = np.append(count, [str('0')  +str(i)])
    for i in range(1000, 10000):
        count = np.append(count, [str(i)])
    return count

def mpi_initial(n_tot, particles, my_rank, p):
    # initialize the particles with Monte-Carlo like method, first try with pure random condition
    #if my_rank == 0:
    #    for source in range(1,p):
    #        add = comm.recv(source=source)
    #        particles = np.append(particles, add)

    #else :
    add = []
    for i in range(int(n_tot/(p-1))):
        check = 1
        while (check):
            position = np.array([rd.uniform(-500, 500), rd.uniform(-500,500), rd.uniform(-500,500)])
            if (position[0]**2+position[1]**2+position[2]**2) <= 250000:
                check = 0
        velocity = np.zeros(3)
        add = np.append(add, [Particle(1e10, position, velocity)])
    #comm.send(add, dest=0)
    #print(str(my_rank)+" finish adding "+str(len(add))+" elements to particle list")
    return add

def initial_solar(n_tot, particles):
    # initialize the particles with Monte-Carlo like method, first try with pure random condition
    # sun as potential, gives 8 planets and n-8 asteroids
    # first generate asteroids
    each_mass = 2.39e24/(n_tot-8)
    for i in range(n_tot-8):
        mass = rd.uniform(0.5*each_mass, 1.5*each_mass)
        r = rd.uniform(2*AU, 3.5*AU)
        theta = rd.uniform(-np.pi, np.pi) # -pi to pi
        phi = rd.uniform(-0.34906585, 0.34906585) # +- 20 degree
        position = np.array([r*math.cos(theta)*math.cos(phi),
        r*math.sin(theta)*math.cos(phi), 
        r*math.sin(phi)])
        v = np.sqrt(G*Msun/r)
        velocity = np.array([-v*math.sin(theta), v*math.cos(theta), 0.0])
        particles = np.append(particles, [Particle(mass, position, velocity)])
    # then add 8 planets (assume at same plane)
    planet_mass = np.array([3.3e26 ,4.87e27 ,5.97e27 ,6.42e26,1.898e30,5.68e29,8.68e28 ,1.02e29])
    planet_r    = np.array([5.79e12,1.082e13,1.496e13,2.28e13,7.78e13,1.432e14,2.867e14,4.515e14])
    planet_vel  = np.array([4.74e6,3.5e6,2.98e6,2.41e6,1.31e6,9.7e5,6.8e5,5.4e5])
    for i in range(8):
        theta = rd.uniform(-np.pi, np.pi)
        mass = planet_mass[i]
        r = planet_r[i]
        v = planet_vel[i]
        position = np.array([r*math.cos(theta),r*math.sin(theta),0.0])
        velocity = np.array([-v*math.sin(theta),v*math.cos(theta),0.0])
#        print(mass)
#        print(position)
#        print(velocity, flush=True)
        particles = np.append(particles, [Particle(mass, position, velocity)])
#        print(particles[i+n_tot-8].pos,flush=True)
#    print(particles.size, flush=True)
    return particles



