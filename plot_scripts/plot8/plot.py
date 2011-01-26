
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

def calc( h0_in, h_in , l ):
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.l = l
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    #h0_in = numpy.array(h0_in)
    #h_in = numpy.array(h_in)
    dmax = 100
    data = []
    data.append(range(1,dmax))
    #
    for h0 in h0_in:
        for h in h_in:
            ising1D.system.h = h
            ising1D.system.h0 = h0
            xxcorrelation1 = []
            xxcorrelation2 = []
            dvect = []
            #
            energy = system.E0(h0,h,l)
            temperature = ising1D.canonical.find_temperature( energy , 0.000001 )
            #
            #print "effective temperature",temperature
            #
            for d in range(1,dmax):
                xxcorrelation1.append( ising1D.canonical.thermal_xxcorrelation(d, temperature))
                xxcorrelation2.append(ising1D.quench.long_time_xxcorrelation( d ) )
                
            data.append( xxcorrelation1 )
            data.append( xxcorrelation2 )
    #
    numpy.savetxt( datadir+'/data.dat', numpy.array(data).transpose() ,fmt='%.18e')
    #
    plot()
    #
    return

def plot():
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir + '/data.dat')
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    #
    for i in range(2,len(data[1,:]),2):
        ax.plot( data[:,0], numpy.log10(numpy.abs(data[:,i])) ,".",ms=4 ,label="long time")
        ax.plot( data[:,0], numpy.log10(numpy.abs(data[:,i-1])) ,".",ms=4 ,label="thermal")
    #
    #
    ax.set_ylim(-3,0)
    ax.legend(loc=0)
    #
    #
    return
