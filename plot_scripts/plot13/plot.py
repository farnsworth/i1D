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

def calc( h0, h , dlist,tmin=0.0,tmax=100.0 , l=512,npoints=1000 ):
    """
    h0          initial h
    h           quench h
    dlist       list of correlation distance
    tmin = 0.0  minimum time
    tmax = 100  maximum time
    l = 512     system length
    npoints = 1000 number of time points
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
    
    tot_line = []
    tot_line2 = []
    tot_line3 = []

    time = numpy.arange(tmin,tmax,(tmax-tmin)/numpy.float(npoints))

    for d in dlist:
        print "d",d
        line = []
        line2 = []
        for t in time:
        #print ising1D.quench.time_xxcorrelation(t, d )
            xx = ising1D.quench.time_xxcorrelation(t, d )
            line.append( numpy.real(xx) )
            line2.append( numpy.imag(xx) )
            if (numpy.abs(numpy.imag(xx))>1.0e-6):
                print "warning"
        tot_line.append(line)
        tot_line2.append(line2)

        tot_line3.append(numpy.ones(numpy.shape(line))*ising1D.quench.long_time_xxcorrelation( d ))
    #
    tot_line.insert(0,time)
    tot_line2.insert(0,time)
    tot_line3.insert(0,time)

    numpy.savetxt( datadir+'/data1.dat', numpy.array(tot_line).transpose(),fmt='%.18e')

    numpy.savetxt( datadir+'/data3.dat', numpy.array(tot_line3).transpose(),fmt='%.18e'  )
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
    data1 = numpy.loadtxt( datadir + '/data1.dat')
    data3 = numpy.loadtxt( datadir + '/data3.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    
    for i in range(1,len(data1[0,:])):
        ax.plot(data1[:,0], numpy.abs(data1[:,i]) ,"-",lw=2)
        #ax.plot(data1[:,0], numpy.abs(data1[:,i]-data3[0,i])/data3[0,i] ,"-",lw=2)
    #    ax.plot(data3[:,0], numpy.abs(data3[:,i]) ,"--",lw=2)
    #
    return

def plot_append():
    #
    ax = pp.gca()
    #
    data1 = numpy.loadtxt( datadir + '/data1.dat')
    data3 = numpy.loadtxt( datadir + '/data3.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    
    for i in range(1,len(data1[0,:])):
        ax.plot(data1[:,0], numpy.abs(data1[:,i]-data3[0,i])/data3[0,i] ,"-",lw=2)
    #    ax.plot(data3[:,0], numpy.abs(data3[:,i]) ,"--",lw=2)
    #

    pp.draw()

    return
