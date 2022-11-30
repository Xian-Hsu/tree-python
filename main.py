import numpy as np
import matplotlib.pyplot as plt
import sys
import os
#sys.path.append(os.path.abspath("/data/xianhsu/others/computational_astrophysics/project"))
from classes import *
from tree_structure import *
from search_tree import *
from visualize import *
from evolve import *
from initial import *
from energy import *
from IO import *
from constants import *

time = 0.0
maxtime = 10000*day
dt = day
step = int(0)
max_step = int(maxtime/dt)
plot_interval = int(1e1)
output_plotdir = np.array(["sliceplot_yz","sliceplot_xy"])
#path = "/data/xianhsu/others/computational_astrophysics/project/plots/"
#file_path = path + output_plotdir
for i in range(2):
    os.makedirs(output_plotdir[i], exist_ok = True)
textfile_name = "data.txt"

n_tot = int(1e3)
n_crit = int(1)
box_length = 70*AU
max_level = 50
particles = []
particles = initial_solar(n_tot, particles)
count = initial_count()
plot_id = int(0)

while (step <= max_step):
    tree = tree_structure(particles, max_level, n_tot, box_length, n_crit)
    for i in range(n_tot):
        rk4(particles, tree, i, time, dt)
    if (step % plot_interval == 0):
        tree_1d = tree_1dlist(tree)
        save_yzsliceplot(tree_1d, box_length, particles, count[plot_id], time) # TODO add time
        save_xysliceplot(tree_1d, box_length, particles, count[plot_id], time) # TODO add time
        E, K, U = cal_energy(particles, tree, n_tot, time)
        L = cal_z_angular_momentum(particles, n_tot)
        write_data(E, K, U, L, plot_id, time, textfile_name)
        plot_id += 1
    time += dt
    step += 1
#tree_1d = tree_1dlist(tree)
#show_yzsliceplot(tree_1d, box_length, particles)
