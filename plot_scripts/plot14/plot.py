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

def calc( h0, h , time ,dmax=100 , l=512 ):
    """
    h0         initial h
    h          quench  h
    time       list of times to perform calculation
    dmax = 100 maximum value of d
    l = 512    system length
    """
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.l = l
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    ising1D.system.h = h
    ising1D.system.h0 = h0
    
    line = []
    line.append(range(2,dmax))

    for t in time:
        line_temp = []
        for d in range(2,dmax):
            xx = ising1D.quench.time_xxcorrelation(t, d )
            line_temp.append( numpy.real(xx) )
            if (numpy.abs(numpy.imag(xx))>1.0e-6):
                print "warning"
        line.append(line_temp)

    line.append([ising1D.quench.long_time_xxcorrelation( d ) for d in line[0] ])        
    #
    numpy.savetxt( datadir+'/data.dat', numpy.array(line).transpose() ,fmt='%.18e')
    #
    plot( time )
    #
    return

def plot( time = [] ):
    """
    It plots the results of calc reading the file data.dat
    """
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir + '/data.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    #
    numdata = len(data[0,:]) - 2
    for i in range(1,numdata+1):
        ax.plot(data[:,0], numpy.real(data[:,i]) ,"-",lw=2,label="t = "+`time[i-1]`)
    #
    ax.plot(data[:,0], numpy.abs(data[:,numdata+1]) ,"--",lw=2,label=r"$t_0 \rightarrow \infty$" )
    #
    ax.legend(loc=0)
    #
    return
