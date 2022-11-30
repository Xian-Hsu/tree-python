import numpy as np

class Particle:
    def __init__(self, mass, position, velocity):
        self.mas = mass
        self.pos = position
        self.vel = velocity
        #self.acc = acceleration

class Cell:
    def __init__(self, position, length): 
        self.mas = 0         # cell mass
        self.pos  = position # minimum value of three axises of cell
        self.apos = np.zeros(3) # weighted average position of cell
        self.len = length    # length of cell
        self.par = []        # index of particles inside cell
        self.sub = []        # subcells of the cell

