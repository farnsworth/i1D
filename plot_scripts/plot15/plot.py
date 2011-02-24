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

def calc( h0, h , d,tmin=0.0,tmax=100.0 , l=512 ):
    """
    h0          initial h
    h           quench h
    d           correlation distance maximum time
    tmax = 100  maximum time
    l = 512     system length
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
    line2 = []

    time = numpy.arange(tmin,tmax,(tmax-tmin)/1000.0)

    for t in time:
#        xx = ising1D.quench.time_biaj(t, d )
        xx = ising1D.quench.time_coeff(t )
        xx2 = ising1D.quench.time_aiaj(t, d )
        line.append( numpy.real(xx) )
        line2.append( numpy.imag(xx2) )
        if (numpy.abs(numpy.imag(xx))>1.0e-6):
            print "warning"

#    line3 = numpy.ones(numpy.shape(line))*ising1D.quench.long_time_coeff( d )
    #
    numpy.savetxt( datadir+'/data.dat', numpy.array([time,line,line2]).transpose() ,fmt='%.18e')
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
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    
    ax.plot(data[:,0], numpy.real(data[:,1]) ,"-",lw=2)
    #ax.plot(data[:,0], numpy.abs(data[:,3]) ,"--",lw=2)
    #
    fig2 = pp.figure()
    ax2  = fig2.add_subplot(111)
    
    ax2.plot(data[:,0], numpy.real(data[:,2]) ,"-",lw=2)
    return

def plot_append():
    
    ax = pp.gca()

    data = numpy.loadtxt( datadir + '/data.dat')

    ax.plot(data[:,0], numpy.real(data[:,1]) ,"-",lw=2)
    ax.plot(data[:,0], numpy.abs(data[:,3]) ,"--",lw=2)

    pp.draw()

    return
