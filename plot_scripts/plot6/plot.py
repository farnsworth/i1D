
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
    npoints = 150
    #
    h0max = 4.0
    h0min = 0.0    
    #
    gsen = ising1D.system.egs()
    emax = max( -h,-1.0 )
    delta = (emax - gsen )/numpy.float(npoints)
    #
    # ... thermal averages
    #
    en = []
    temp = []
    xxcorrelation = []

    for i in range(1,npoints):
        en.append( gsen + delta*numpy.float(i) )
        temp.append( ising1D.canonical.find_temperature( en[i-1], 0.000001 ) )
        xxcorrelation.append( ising1D.canonical.thermal_xxcorrelation(d, temp[i-1]))
    
    data = numpy.column_stack( ( en, temp, xxcorrelation ) )
    numpy.savetxt( datadir+'/data_thermal.dat', data ,fmt='%.18e')
    #
    # ... time average
    #
    en2 = []
    xxcorrelation2 = []
    delta = (h0max - h0min)/numpy.float(npoints)
    #
    for i in range(0,npoints):
        h0 = delta*i + h0min
        ising1D.system.h0 = h0
        en2.append( system.E0(h0,h,l) )
        xxcorrelation2.append( ising1D.quench.long_time_xxcorrelation( d ) )
    #
    data = numpy.column_stack( ( en2, xxcorrelation2 ) )
    numpy.savetxt( datadir+'/data_time.dat', data ,fmt='%.18e')
    #
    plot()
    #
    return

def plot():
    #
    print "... Doing the plot ..."
    #
    data_thermal = numpy.loadtxt( datadir + '/data_thermal.dat')
    data_time = numpy.loadtxt( datadir + '/data_time.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.plot(data_thermal[:,0], data_thermal[:,2] ,".",ms=4, label="thermal" )
    ax.plot(data_time[:,0],data_time[:,1],".",ms=4,label="time")
    #
    return
