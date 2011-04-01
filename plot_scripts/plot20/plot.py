
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

def calc( h0, h, deltae , l =30,filename="data.dat" ):
    """Do a plot of microcanonical fluctuation and domain of xx coorelation
    as a function of l.
    
    Usage:
          calc(d, h0, h, deltae, xxbin, xxmax, l)

    Where:
          d         --> correlation distance
          h0        --> initial or transverse magnetization
          h         --> transverse magnetization
          deltae    --> microcanonical energy window
          xxbin=100 --> number of xx bins
          xxmax=1.0 --> maximum value of xxcorrelation
          l=30      --> size of the system
    """
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.h = h
    ising1D.system.h0 = h0
    ising1D.system.l = l
    ising1D.system.datadir[:len(datadir)] = datadir
    xxmax = 1.0
    xxbin = 100
    energy = system.E0(h0,h,l)
    #
    data = []
    #
    # ... doing calculation
    for d in range(2,30,1):
        sigma = ising1D.microcanonical.xxcorrelation_sigma_mc(d,energy, deltae, xxmax, xxbin)
        data.append( [d,ising1D.microcanonical.obsmax.copy(),ising1D.microcanonical.obsmin.copy(),sigma ] )
        print "max",ising1D.microcanonical.obsmax
    #
    # ... saving data
    #
    numpy.savetxt( datadir+'/'+filename, data, fmt='%.18e')
    #
    plot(filename)
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
    ax.plot( data[:,0],data[:,1],label=filename)
    ax.plot( data[:,0],data[:,2],label=filename)
    #ax.plot( data[:,0],data[:,3])
    #
    return

#
#def append_plot(ax,filename):
#    data = numpy.loadtxt( datadir+'/'+filename)
#    ax.plot( data[:,0],data[:,3],label=filename)
#
#def plot1(h0,h):
#    calc( 4, h0, h, 0.05, lmax=60, filename="data005.dat")
#    calc( 4, h0, h, 0.1, lmax=60, filename="data01.dat")
#    calc( 4, h0, h, 0.2, lmax=60, filename="data02.dat")
#    fig = pp.figure()
#    ax  = fig.add_subplot(111)
#    append_plot(ax, "data005.dat")
#    append_plot(ax, "data01.dat")
#    append_plot(ax, "data02.dat")
#    ax.legend()
