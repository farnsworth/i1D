
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

def calc( d, h0, h, deltae, l=30,filename="data.dat" ):
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.h0 = h0
    ising1D.system.h = h
    ising1D.system.l = l
    energy = system.E0(h0,h,l)
    #
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    random_seed = numpy.random.random_integers(2024, size=8)
    #
    xxmax = 1.0
    xxmin = -1.0
    #
    ising1D.wl_rho.accuracy = 1.0e-3
    #
    out = ising1D.wl_rho.wl_rhodos( d, energy, deltae, xxmin,xxmax,500,random_seed,False)
    #
    print "exact calculation"
    ising1D.microcanonical.xxcorrelation_sigma_mc(d, energy, deltae, xxmax, 500)
    #
    if out not in (0,2):
        return
    #
    # ... saving data
    #
    #numpy.savetxt( datadir+'/'+filename, ising1D.wl_rho.logdosf.transpose(), fmt='%.18e')
    numpy.savetxt( datadir+'/'+filename, ising1D.wl_rho.logdosf, fmt='%.18e')
    numpy.savetxt( datadir+'/'+"data_exact.dat", ising1D.microcanonical.dist, fmt='%.18e')
    #
    plot()
    #
    return

def plot(filename="data.dat"):
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir+'/'+filename)
    data_exact = numpy.loadtxt( datadir+'/'+"data_exact.dat")
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.plot(data[:,0],data[:,1],'o',label='numerical' )
    for i in range(len(data_exact[:,0])):
        if (data_exact[i,1]>0):
            alpha = numpy.log(data_exact[i,1])
            break

    ax.plot(data_exact[:,0],numpy.log(data_exact[:,1])-alpha,'o',label='exact')
    ax.legend(loc=0)
    #
    return
