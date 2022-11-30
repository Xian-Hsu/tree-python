import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from datetime import datetime
#import os
#import sys
#sys.path.append(os.path.abspath("/data/xianhsu/others/computational_astrophysics/project"))
from classes import *
from tree_structure import *
from search_tree import *
from constants import *

def show_3dtree(tree_1d):
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    xdata = []
    ydata = []
    zdata = []
    for i in range(len(tree_1d)):
        xdata = np.append(xdata, tree_1d[i].pos[0])
        ydata = np.append(ydata, tree_1d[i].pos[1])
        zdata = np.append(zdata, tree_1d[i].pos[2])
        #print('xdata, ydata, zdata is %.2f, %.2f, and %.2f' % (xdata[i], ydata[i], zdata[i]))
        #print(i)
    ax.scatter3D(xdata, ydata, zdata)#, c=zdata, cmap='Greens')
    #plt.show()
    plt.clf()

def show_yzsliceplot(tree_1d, box_length, particles):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(len(tree_1d)):
        if tree_1d[i].pos[0] == 0:
            #print(tree_1d[i].pos)
            ax.plot([tree_1d[i].pos[1], tree_1d[i].pos[1]+tree_1d[i].len],
            [tree_1d[i].pos[2], tree_1d[i].pos[2]], 'k')
            ax.plot([tree_1d[i].pos[1], tree_1d[i].pos[1]],
            [tree_1d[i].pos[2], tree_1d[i].pos[2]+tree_1d[i].len], 'k')
            ax.plot([tree_1d[i].pos[1], tree_1d[i].pos[1]+tree_1d[i].len],
            [tree_1d[i].pos[2]+tree_1d[i].len, tree_1d[i].pos[2]+tree_1d[i].len], 'k')
            ax.plot([tree_1d[i].pos[1]+tree_1d[i].len, tree_1d[i].pos[1]+tree_1d[i].len],
            [tree_1d[i].pos[2], tree_1d[i].pos[2]+tree_1d[i].len], 'k')
            ax.plot(tree_1d[i].apos[1], tree_1d[i].apos[2], 'bx')
            if len(tree_1d[i].par) != 0:
                ax.plot(particles[int(tree_1d[i].par[0])].pos[1], particles[int(tree_1d[i].par[0])].pos[2], 'r.')
    ax.set_xlim(-box_length/2, box_length/2)
    ax.set_ylim(-box_length/2, box_length/2)
    #plt.show()
    plt.clf()

def save_yzsliceplot(tree_1d, box_length, particles, number, time):
    fig = plt.figure()
    plt.clf()
    ax = fig.add_subplot(111)
    for i in range(len(tree_1d)):
        if tree_1d[i].pos[0] == 0:
            #print(tree_1d[i].pos)
            ax.plot([tree_1d[i].pos[1], tree_1d[i].pos[1]+tree_1d[i].len],
            [tree_1d[i].pos[2], tree_1d[i].pos[2]], 'k')
            ax.plot([tree_1d[i].pos[1], tree_1d[i].pos[1]],
            [tree_1d[i].pos[2], tree_1d[i].pos[2]+tree_1d[i].len], 'k')
            ax.plot([tree_1d[i].pos[1], tree_1d[i].pos[1]+tree_1d[i].len],
            [tree_1d[i].pos[2]+tree_1d[i].len, tree_1d[i].pos[2]+tree_1d[i].len], 'k')
            ax.plot([tree_1d[i].pos[1]+tree_1d[i].len, tree_1d[i].pos[1]+tree_1d[i].len],
            [tree_1d[i].pos[2], tree_1d[i].pos[2]+tree_1d[i].len], 'k')
            ax.plot(tree_1d[i].apos[1], tree_1d[i].apos[2], 'bx')
            if len(tree_1d[i].par) != 0:
                ax.plot(particles[int(tree_1d[i].par[0])].pos[1], particles[int(tree_1d[i].par[0])].pos[2], 'r.')
    ax.set_xlim(-box_length/2, box_length/2)
    ax.set_ylim(-box_length/2, box_length/2)
    ax.set_title("time = %i days" %int(time/day))
    fig.savefig("sliceplot_yz/"+number+".png")
    #plt.figure().clear()
    plt.close(fig)
    print("complete figure number " + number + " at " + str(datetime.now()), flush=True)
    #plt.cla()
    #plt.clf()

def save_xysliceplot(tree_1d, box_length, particles, number, time):
    fig = plt.figure()
    plt.clf()
    ax = fig.add_subplot(111)
    for i in range(len(tree_1d)):
        if tree_1d[i].pos[2] == 0:
            #print(tree_1d[i].pos)
            ax.plot([tree_1d[i].pos[0], tree_1d[i].pos[0]+tree_1d[i].len],
            [tree_1d[i].pos[1], tree_1d[i].pos[1]], 'k')
            ax.plot([tree_1d[i].pos[0], tree_1d[i].pos[0]],
            [tree_1d[i].pos[1], tree_1d[i].pos[1]+tree_1d[i].len], 'k')
            ax.plot([tree_1d[i].pos[0], tree_1d[i].pos[0]+tree_1d[i].len],
            [tree_1d[i].pos[1]+tree_1d[i].len, tree_1d[i].pos[1]+tree_1d[i].len], 'k')
            ax.plot([tree_1d[i].pos[0]+tree_1d[i].len, tree_1d[i].pos[0]+tree_1d[i].len],
            [tree_1d[i].pos[1], tree_1d[i].pos[1]+tree_1d[i].len], 'k')
            ax.plot(tree_1d[i].apos[0], tree_1d[i].apos[1], 'bx')
            if len(tree_1d[i].par) != 0:
                ax.plot(particles[int(tree_1d[i].par[0])].pos[0], particles[int(tree_1d[i].par[0])].pos[1], 'r.')
    ax.set_xlim(-box_length/2, box_length/2)
    ax.set_ylim(-box_length/2, box_length/2)
    ax.set_title("time = %i days" %int(time/day))
    fig.savefig("sliceplot_xy/"+number+".png")
    #plt.figure().clear()
    plt.close(fig)
    print("complete figure number " + number + " at " + str(datetime.now()), flush=True)
    #plt.cla()
    #plt.clf()

