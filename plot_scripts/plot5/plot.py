
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
from  time import sleep,time

from matplotlib import rc
rc('text', usetex=True)

def calc( d, h, l ):
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.l = l
    ising1D.system.h = h
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    # ... doing calculation in all the space
    print " ...All space calculation..."
    ising1D.exact.corralpha_calc(d)
    # ... saving data
    numpy.savetxt( datadir+'/data.dat', ising1D.exact.real_array.transpose() , fmt='%.18e')
    #
    # ... doing calculation in the sub space with symmetry k,-k
    print " ...Reduced space calculation..."
    ising1D.exact.corralpha_calc_red(d)
    # ... saving data
    numpy.savetxt( datadir+'/data_red.dat', ising1D.exact.real_array.transpose() , fmt='%.18e')
    #
    plot()
    #
    return

def plot():
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir+'/data.dat')
    data_red = numpy.loadtxt( datadir+'/data_red.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.plot(data[:,0], data[:,1] ,".",ms=4, label="numerical" )
    ax.plot(data_red[:,0], numpy.abs(data_red[:,1]) ,"x",ms=6, label="numerical" )

    return
