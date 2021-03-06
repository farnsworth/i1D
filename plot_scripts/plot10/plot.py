"""
It does a plot of xx-correlation for an absent quench h=h_0
"""
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
#from  time import sleep,time

from matplotlib import rc
rc('text', usetex=True)

def calc( h_vect, d_max ,l ):
    """ 
    Perform calculation of xx-correlation for an absent quench.
    Inputs:
       - h_vect : list with values of h
       - d_max : maximum distance
       - l : dimension of the system
    """
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.l = l
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    #
    xxcorrelations = []
    xxcorrelations.append(range(1,d_max+1))
    #
    #
    xxcorrelations2 = []
    #
    xxcorrelations2.append(range(1,d_max+1))
    #
    #
    d_max_fit = d_max
    d_min_fit = 80
    #
    #
    for h in h_vect:
        
        h_limit = 0.0

        if (h < 1.0):
            h_limit = (1.0-h*h)**(0.25)
        
        ising1D.system.h = h
        ising1D.system.h0 = h
        line = []
        line2 = []
        for d in range(1,d_max+1):
            line.append( ising1D.quench.long_time_xxcorrelation( d ) - h_limit )
            line2.append( ising1D.canonical.thermal_xxcorrelation(d, 0.0 ))
 
        xxcorrelations.append(line)
        xxcorrelations2.append(line2)

        data = numpy.log(line[d_min_fit-1:d_max_fit])
        xdata = numpy.arange(d_min_fit,d_max_fit+1)

#        print len(data)
#        print len(xdata)

        p = numpy.polyfit(xdata,data,1)

        print "h",h,"length, theoretical",-numpy.log(numpy.log(h)),"from fit",-numpy.log(-p[0])

    #
    numpy.savetxt( datadir+'/data.dat', numpy.array(xxcorrelations).transpose() ,fmt='%.18e')
    numpy.savetxt( datadir+'/data2.dat', numpy.array(xxcorrelations2).transpose() ,fmt='%.18e')
    #
    plot()
    #
    return

def plot():
    """
    It plots the results of calc reading the file data.dat
    """
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir + '/data.dat')
    data2 = numpy.loadtxt( datadir + '/data2.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    for i in range(len(data[1,:])-1):
        ax.plot(data[:,0], numpy.log(data[:,i+1]) ,".",ms=8 )
        ax.plot(data2[:,0], numpy.log(data2[:,i+1]) ,".",ms=3 )
    #
    return
