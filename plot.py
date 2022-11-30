import matplotlib.pyplot as plt
import numpy as np

# E K U L T

colors = np.array(['b','k','r','m'])
file_name = "data.txt"

with open(file_name) as textFile:
    lines = [line.split() for line in textFile]
data = np.array(lines[1:], dtype = float)
E = np.array(data[:,0], dtype = float)
K = np.array(data[:,1], dtype = float)
U = np.array(data[:,2], dtype = float)
L = np.array(data[:,3], dtype = float)
T = np.array(data[:,4], dtype = float)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.semilogy(T, abs(E), colors[0])
ax.semilogy(T,     K , colors[1])
ax.semilogy(T, abs(U), colors[2])
ax.semilogy(T, abs(L), colors[3])
ax.legend(['abs(total energy)','kinetic energy','abs(potential energy)','abs(z angular momentum)'])
fig.savefig("energy.png")



