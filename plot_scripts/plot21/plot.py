
import sys
import os

datadir = os.path.dirname(os.path.abspath(__file__))

sys.path.append(datadir+'/../../libs/')
sys.path.append(datadir+'/../../')

import numpy
import math
import ising1D
import utilities
import system
import matplotlib.pyplot as pp
import matplotlib.cm as cm
from  time import sleep,time

from matplotlib import rc
rc('text', usetex=True)

def calc( h, l=1000,nstates=100,method="random",filename="data.dat" ):
    """Do a plot of microcanonical fluctuation and domain of xx coorelation
    as a function of l.
    
    Usage:
          calc(d, h0, h, deltae, xxbin, xxmax, l)

    Where:
          h                   --> transverse field
          l                   --> length of the system
          nstates=100         --> number of states 
          method="random"     --> method: "random" or "sequential"
          filename="data.dat" --> output filename
    """
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    state = numpy.random.rand(l/2)<0.5
    #
    #
    #energy = system.E0(h0,h,l)
    #
    data = []
    datatot = []
    #
    fig = pp.figure()
    fig_sign = pp.figure()
    ax  = fig.add_subplot(111)
    ax_sign = fig_sign.add_subplot(111) 
    egs = ising1D.system.egs()
    #
    # ... doing calculation
    #
    for i in range(0,nstates):
        if (method=="random"):
            state = numpy.random.rand(l/2)<0.5
        else:
            state = numpy.zeros(l/2,dtype=bool)
            binnum = bin(i)
            for index in range( len(binnum)-1, 1, -1 ):
                if (binnum[index]=='1'):
                    state[len(binnum)-1-index]=True

        data = []
        for d in range(2,100,1):
            val = ising1D.system.xxcorrelation_red(d,state,l/2)
            data.append( [d,val ] )

        arr = numpy.array(data)
        ax.plot(arr[:,0],numpy.log( numpy.abs(arr[:,1])))
        ax_sign.plot(arr[:,0],numpy.sign( arr[:,1] ) )
        #nminfit = 30
        #nmaxfit = 100
        #rate = m(arr[nminfit:nmaxfit,0],numpy.log(numpy.abs(arr[nminfit:nmaxfit,1])) )
        #energy = egs + ising1D.system.state_ex_energy(state,l/2)
        #datatot.append( [energy,-1.0/rate] )
    #
    ax_sign.set_ylim(-1.2,1.2)
    #
    # ... saving data
    #
    #numpy.savetxt( datadir+'/'+filename, datatot, fmt='%.18e')
    #
    #plot(filename)
    #
    return

def plot(filename="data.dat"):
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir+'/'+filename)
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    #
    ax.plot( data[:,0],data[:,1],"o",label=filename)
    #
    return


def m(xdata,ydata):
    sx = xdata.sum()
    sy = ydata.sum()
    sxx = (xdata*xdata).sum()
    sxy = (xdata*ydata).sum()
    n = len(xdata)
    return (n*sxy-sx*sy)/(n*sxx-sx*sx)
