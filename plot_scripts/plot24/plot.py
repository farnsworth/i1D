
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

def calc( d, h0, h, deltae,filename="data" ):
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.h0 = h0
    ising1D.system.h = h
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
    ising1D.wl_rho.accuracy = 1.0e-5
    data = []
    #
    energy = system.E0(h0,h,100)
    #
    for i in range(20,130,4):
        print "\nL = ",i
        ising1D.system.l = i
        out = ising1D.wl_rho.wl_rhodos( d, energy, deltae, xxmin,xxmax,nbin,random_seed,False)
        #
        #
        if out not in (0,2):
            print "...Simulation stopped..."
            continue

        res = utilities.HistoFromLog(ising1D.wl_rho.logdosf[:,1],delta)
        sigma = calc_sigma( ising1D.wl_rho.logdosf[:,0],res, delta )
        print "sigma",sigma
        data.append([i,sigma,out])
    #
    # ... saving data
    #
    numpy.savetxt( datadir+'/'+filename+'.dat', data, fmt='%.18e')
    #
    plot()
    #
    return

def plot(filename="data.dat"):
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir+'/'+filename)
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.plot(data[:,0],data[:,1],'o')
    #
    return

def calc_sigma(xdata, ydata, delta):
    tot = mean_x = mean_x2 = 0.0
    for k in range(len( xdata)):
        mean_x += xdata[k]*ydata[k]*delta
        mean_x2 += xdata[k]*xdata[k]*ydata[k]*delta
        tot += delta*ydata[k]

    return numpy.sqrt( mean_x2 - mean_x*mean_x )



#### generate plot ####
def gen_plot(d):
    h=0.5
    h0=0.9

    deltae = 0.1
    filename="data-01"
    calc(d,h0,h,deltae,filename)

    deltae = 0.05
    filename="data-005"
    calc(d,h0,h,deltae,filename)

    # deltae too small for this order of L
    #deltae = 0.01
    #filename="data-001"
    #calc(d,h0,h,deltae,filename)

def plot2():
    data1 = numpy.loadtxt( datadir+'/data-01.dat')
    data2 = numpy.loadtxt( datadir+'/data-005.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.plot(data1[:,0],data1[:,1],'o',label = '$\delta = 0.1$')
    ax.plot(data2[:,0],data2[:,1],'o',label = '$\delta = 0.05$')
    ax.legend(loc=0)
