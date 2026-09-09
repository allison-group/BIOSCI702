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
    y = []
    for xVal in x:
        sigOverR6 = np.power(ensemble.sigma/xVal, 6)
        sigOverR12 = sigOverR6*sigOverR6        
        y.append( 4 * ensemble.epsilon * (sigOverR12 - sigOverR6))
    df = pd.DataFrame()
    df['distance'] = x
    df['LJ potential energy'] = y
    df.plot(x = 'distance', y = 'LJ potential energy', title = 'Lennard-Jones potential', ylabel = 'potential energy')
    #plt.show()

def plot2LJpotentials(ensemble1, ensemble2, ptCount=200):
    # first handle ensemble 1
    xlowerlim1 = ensemble1.sigma - 0.05
    x1 = np.linspace(xlowerlim1, ensemble1.cutoff, ptCount)
    y1 = []
    for xVal in x1:
        sigOverR6 = np.power(ensemble1.sigma/xVal, 6)
        sigOverR12 = sigOverR6*sigOverR6        
        y1.append( 4 * ensemble1.epsilon * (sigOverR12 - sigOverR6))
    df1 = pd.DataFrame()
    df1['distance'] = x1
    df1['LJ1'] = y1
    # then handle ensemble 2
    xlowerlim2 = ensemble2.sigma - 0.05
    x2 = np.linspace(xlowerlim2, ensemble2.cutoff, ptCount)
    y2 = []
    for xVal in x2:
        sigOverR6 = np.power(ensemble2.sigma/xVal, 6)
        sigOverR12 = sigOverR6*sigOverR6        
        y2.append( 4 * ensemble2.epsilon * (sigOverR12 - sigOverR6))
    df2 = pd.DataFrame()
    df2['distance'] = x2
    df2['LJ2'] = y2
    # now combine them
    df = pd.concat([df1.LJ1, df2.LJ2], axis=1)
    dfForPlot = pd.concat([df1.distance, df], axis=1)
    dfForPlot.columns = ['distance','LJ1','LJ2']
    dfForPlot.plot(x='distance', y=['LJ1','LJ2'], title='Lennard-Jones potentials', ylabel = 'potential energy')
    plt.show()

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
    
    
