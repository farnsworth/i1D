
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

def calc( h0, h, deltae, l=200, filename="data" ):
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.h0 = h0
    ising1D.system.h = h
    ising1D.system.l = l
    #
    ising1D.system.datadir[:len(datadir)+1] = datadir+'/'
    #
    random_seed = numpy.random.random_integers(2024, size=8)
    #
    xxmax = 1.0
    xxmin = -1.0
    #
    nbin = 1000
    delta = (xxmax-xxmin)/float(nbin)
    #
    ising1D.wl_rho.accuracy = 1.0e-2
    data = []
    #
    energy = system.E0(h0,h,l)
    #
    for i in range(2,50):
        #
        out = ising1D.wl_rho.wl_rhodos( i, energy, deltae, xxmin,xxmax,nbin,random_seed,False)
        #
        if (out != 0):
            print "...Simulation stopped..."
            return

        data.append( [ i, ising1D.wl_rho.rhomin_out, ising1D.wl_rho.rhomax_out ] )
    #
    # ... saving data
    #
    numpy.savetxt( datadir+'/'+filename+'.dat', data, fmt='%.18e')
    #
    plot()
    #
    return


def plot(filename="data"):
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir+'/'+filename+'.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    #
    ax.plot(data[:,0],data[:,1],'o',label='max')
    ax.plot(data[:,0],data[:,2],'o',label='min')
    #
    return
