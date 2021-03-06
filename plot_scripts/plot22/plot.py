
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

def calc( d, h, l=30,filename="data.dat" ):
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    #
    gsen = ising1D.system.egs()
    #
    random_seed = numpy.random.random_integers(2024, size=8)
    #
    emin = gsen
    emax = -gsen
    xxmax = 1.0
    xxmin = -1.0
    #
    ising1D.wl_rho.accuracy = 1.0e-3
    #
    out = ising1D.wl_rho.wl_rhojdos(d,emin,emax,xxmin,xxmax,500,random_seed)
    #
    if out not in (0,2):
        return
    #
    # ... saving data
    #
    #numpy.savetxt( datadir+'/'+filename, ising1D.wl_rho.logdosf.transpose(), fmt='%.18e')
    numpy.savetxt( datadir+'/'+filename, ising1D.wl_rho.domain.transpose(), fmt='%.18e')
    #
    plot([emin,emax],[xxmax,xxmin],filename)
    #
    return

def plot(xlim,ylim,filename="data.dat"):
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir+'/'+filename)
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    cax = ax.imshow( 1-data,origin="lower",interpolation="nearest",
                     extent=(xlim[0],xlim[1],ylim[0],ylim[1]),
                     aspect ='auto',cmap=cm.gray_r)
    cbar = fig.colorbar( cax )
    #cbar.set_label(r"$\log ( n(E,\rho^{xx}_d) + 1 ) $")
    #                                                                         
    #ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    #
    return
