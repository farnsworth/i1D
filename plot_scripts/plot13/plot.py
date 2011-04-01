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

################ 1 --- specific calculation ##################

def generate_plot1():
    calc(2.0,1.25,[20,40,80],tmin=0.0, tmax = 150.0, l=500, npoints=2000 )
    plot_with_zoom1(100,125,[20,40,80] )
    return

def plot_with_zoom1(fr,to, labels):
    
    h = 2.0
    l = 500

    data1 = numpy.loadtxt( datadir + '/data1.dat')
    
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\rho^{xx}_l$")
    ax2 = fig.add_axes( (0.25,0.4,0.35,0.4) )
    ax2.set_xlabel(r"$t$")
    ax2.set_ylabel(r"$\rho^{xx}_l$")

    for i in range(1,len(data1[0,:])):
        line = ax.plot(data1[:,0], numpy.abs(data1[:,i]) ,"-",lw=2,label="l = "+str(labels[i-1]) )[0]
        ax2.plot(data1[:,0], numpy.abs(data1[:,i]) ,"-",lw=2,color=line.get_color())
        ax2.axvline(x=tstar(labels[i-1],l,h),lw=1.5,color=line.get_color(),linestyle="--")
        ax.axvline(x=tstar(labels[i-1],l,h),lw=1.5,color=line.get_color(),linestyle="--")

    ax2.set_xlim(fr,to)
    ax2.set_ylim(top=0.03)
    ax.legend()
    pp.draw()

    return

#########################################################

################ 2 --- specific calculation #############

def generate_plot2():
    calc(0.9,0.5,[20,40,80],tmin=0.0,tmax=300.0,l=500,npoints=300)
    plot_with_zoom2(200,260,[20,40,80] )
    return

def plot_with_zoom2(fr,to, labels):
    
    h = 0.5
    l = 500

    data1 = numpy.loadtxt( datadir + '/data1.dat')
    
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\rho^{xx}_l$")
    ax2 = fig.add_axes( (0.25,0.4,0.35,0.4) )
    ax2.set_xlabel(r"$t$")
    ax2.set_ylabel(r"$\rho^{xx}_l$")

    for i in range(1,len(data1[0,:])):
        line = ax.plot(data1[:,0], numpy.abs(data1[:,i]) ,"-",lw=2,label="l = "+str(labels[i-1]) )[0]
        ax2.plot(data1[:,0], numpy.abs(data1[:,i]) ,"-",lw=2,color=line.get_color())
        ax2.axvline(x=tstar(labels[i-1],l,h),lw=1.5,color=line.get_color(),linestyle="--")
        ax.axvline(x=tstar(labels[i-1],l,h),lw=1.5,color=line.get_color(),linestyle="--")

    ax2.set_xlim(fr,to)
    ax2.set_ylim(top=0.25)
    ax.legend()
    pp.draw()

    return

#########################################################

def tstar(d,l,h):
    if (h>=1):
        return float(l-d)/4.0
    else:
        return float(l-d)/(4.0*h)


################# -- plot of fluctuations ##############

def generate_plot3():
    calc(2.0,1.25,[80,40,20],tmin=0.0, tmax = 1000.0, l=5000, npoints=5000 )
    plot3( [80,40,20] )
    return

def plot3 ( labels ):
    #
    h0 = 2.0
    h = 1.25
    l = 5000
    #
    data1 = numpy.loadtxt( datadir + '/data1.dat')
    data3 = numpy.loadtxt( datadir + '/data3.dat')
    
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\frac{\rho^{xx}_l(t) - <\rho^{xx}_l>_\mathrm{D}}{<\rho^{xx}_l>_\mathrm{D}}$")

    for i in range(1,len(data1[0,:])):
        ax.plot(data1[:,0], numpy.abs((data1[:,i]-data3[:,i])/data3[:,i]) ,"-",lw=2,label="l = "+str(labels[i-1]) )[0]

    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    pp.draw()

    return
