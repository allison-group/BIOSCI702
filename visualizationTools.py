#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 28 14:43:33 2021

@author: Gianmarc Grazioli

"""

import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
    
def makeTrajMovie2D(traj, sideLen, filename = 'LJtraj.gif'):
    fig = plt.figure(figsize=(5,5))
    plt.xlim(0 - .5, sideLen+.5)
    plt.ylim(0 - .5, sideLen+.5)
    graph, = plt.plot([], [], 'o')
    def animate(i):
        x = traj[i][:,0]
        y = traj[i][:,1]
        graph.set_data(x, y)
        return(graph,)
    ani = FuncAnimation(fig, animate, frames=len(traj), blit=True)
    writergif = PillowWriter(fps=30) 
    ani.save(filename, writer=writergif)
    
def makeTrajMovie2DColored(traj, sideLen, filename = 'LJtraj.gif'):
    fig = plt.figure(figsize=(5,5))
    plt.xlim(0 - .5, sideLen+.5)
    plt.ylim(0 - .5, sideLen+.5)
    sphere = 10
    colors = np.random.rand(sideLen)
    graph, = plt.scatter([], [], s=sphere, c=colors)
    def animate(i):
        x = traj[i][:,0]
        y = traj[i][:,1]
        graph.set_data(x, y)
        return(graph,)
    ani = FuncAnimation(fig, animate, frames=len(traj), blit=True)
    writergif = PillowWriter(fps=30)
    ani.save(filename, writer=writergif)

def plotKEtotals(filepath):
    KEdf = pd.read_csv(filepath, header=0)
    KEdf.plot(x='time', y='kinetic energy')
    
def plotPEtotals(filepath):
    PEdf = pd.read_csv(filepath, header=0)
    PEdf.plot(x='time', y='potential energy')
    
def plotKEandPE(filepathKE, filepathPE):
    KEdf = pd.read_csv(filepathKE, header=0)
    PEdf = pd.read_csv(filepathPE, header=0)
    df = pd.concat([KEdf, PEdf.PE], axis=1)
    df.plot(x='time', y=['kinetic energy','potential energy'])
    
def plotTotalEnergy(filepathKE, filepathPE):
    KEdf = pd.read_csv(filepathKE, header=0)
    PEdf = pd.read_csv(filepathPE, header=0)
    df = pd.concat([KEdf.KE, PEdf.PE], axis=1)
    dfForPlot = pd.concat([KEdf.t, df.sum(axis=1)], axis=1)
    dfForPlot.columns = ['time', 'total energy']
    dfForPlot.plot(x='time', y='total energy')
    
def plotKEandPEandTotal(filepathKE, filepathPE): 
    KEdf = pd.read_csv(filepathKE, header=0)
    PEdf = pd.read_csv(filepathPE, header=0)
    df = pd.concat([KEdf.KE, PEdf.PE], axis=1)
    dfForPlot = pd.concat([KEdf.t, df, df.sum(axis=1)], axis=1)
    dfForPlot.columns = ['time', 'kinetic energy', 'potential energy', 'total energy']
    dfForPlot.plot(x='time', y=['kinetic energy','potential energy', 'total energy'])
    
def plotLJpotential(ensemble, ptCount=200):
    #x = np.linspace(ensemble.sigma, ensemble.cutoff, ptCount)
    xlowerlim = ensemble.sigma - 0.05
    x = np.linspace(xlowerlim, ensemble.cutoff, ptCount)
    #print(x[0])
    y = []
    for xVal in x:
        sigOverR6 = np.power(ensemble.sigma/xVal, 6)
        sigOverR12 = sigOverR6*sigOverR6        
        y.append( 4 * ensemble.epsilon * (sigOverR12 - sigOverR6))
    df = pd.DataFrame()
    df['distance'] = x
    df['potential energy'] = y
    df.plot(x = 'distance', y = 'potential energy', title = 'Lennard-Jones potential', ylabel = 'potential energy')
    #plt.show()

def writeXYZ(ensemble, filename, atomName='Ar'):
    f = open(filename, 'w')
    molNum = str(ensemble.nMol) 
    is2D = (len(ensemble.trajectory[0][0]) == 2)
    for step in ensemble.trajectory:
        f.write(molNum+'\n\n')
        for mol in step:
            f.write(atomName+'\t')
            for dim in mol:
                f.write('%.4f    ' % dim)
            if (is2D): 
                f.write("0.0000")
            f.write('\n')
    f.close()     
    
    
