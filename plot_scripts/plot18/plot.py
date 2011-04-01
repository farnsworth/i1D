
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

def calc( d, h, ebin=100, xxbin=100 ,xxmax = 1.0 , l =30 ):
    """Do a plot of 2D density of states in the plane 
       xx-energy (eigenstates with simmetry k,-k).
    
    Usage:
          calc(d, h, xxmax)

    Where:
          d         --> correlation distance
          h         --> transverse magnetization
          ebin=100  --> number of energy bins
          xxbin=100 --> number of xx bins
          xxmax=1.0 --> maximum value of xxcorrelation
          l=30      --> size of the system
    """
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.l = l
    ising1D.system.h = h
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    # ... doing calculation
    ising1D.exact.xxcorralpha_array_calc(d,ebin,xxbin,xxmax)
    #
    # ... saving data
    #
    numpy.savetxt( datadir+'/data.dat', ising1D.exact.array.transpose(), fmt='%.18e')
    #
    plot(ising1D.exact.emax_plot,xxmax,[ising1D.exact.obsmin_plot,ising1D.exact.obsmax_plot])
    #
    return

def plot(emax, xxmax, xxlim ):
    #
    print "... Doing the plot ..."
    #
    data_tmp = numpy.loadtxt( datadir+'/data.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    data = numpy.log(data_tmp + 1.0)
    cax = ax.imshow( data,origin="lower",interpolation="nearest",
                     extent=(-emax,emax,-xxmax,xxmax),
                     aspect ='auto',cmap=cm.gray_r)
    cbar = fig.colorbar( cax )
    cbar.set_label(r"$\log ( n(E,\rho^{xx}_d) + 1 ) $")
    #
    ax.set_xlim(-emax,emax)
    ax.set_ylim( xxlim)
    ax.set_xlabel(r"$\epsilon_\alpha$")
    ax.set_ylabel(r"$\rho^{xx}_d$")
    #
    return
