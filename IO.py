import numpy as np
from classes import *

def write_data(energy, K, U, angular_momentum, plot_id, time, textfile_name):
    # write simulation data such as energy and time into .txt file
    f = open(textfile_name, "a")
    if plot_id == 0:
        f.write("energy kinetic_energy potential_energy angular_momentum time\n")
    f.write("%.5e " %energy)
    f.write("%.5e " %K)
    f.write("%.5e " %U)
    f.write("%.5e " %angular_momentum)
    f.write("%.5e\n" %time)
    f.close()

#def output_particles(particles):
#   # output particles as np.array form


