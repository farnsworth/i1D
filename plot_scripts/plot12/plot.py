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

def exact(h,temp):
    c = 2.0*numpy.sqrt(h)
    m = (1.0-h)/(2.0*h)
    delta = 2.0*numpy.abs(1-h)
    res = numpy.exp(delta/temp)*numpy.sqrt(numpy.pi*c/(2.0*delta*temp))
    return res

def exact2(h,temp):
    c = 2.0*numpy.sqrt(h)
    m = (1.0-h)/(2.0*h)
    delta = 2.0*numpy.abs(1-h)
    xi_old = exact(h,temp)
    res = c*xi_old/(delta*xi_old+c)
    print numpy.log10(c/delta)
    return res

    

def calc( h_vect, temp_vect , d_max ,l ):
    """ 
    Perform calculation of xx-correlation for different temperatures.
    Inputs:
       - h_vect : list with values of h
       - temp_vect : temperature list
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
    xi = []
    xi_theory = []
    #
    #
    d_max_fit = d_max
    d_min_fit = 50
    #
    fig = pp.figure()
    axes = fig.add_subplot(111)
    #
    #
    for h in h_vect:
        #
        #
        #if (h > 1.0):
        #    continue
        
        for temp in temp_vect:
            #
            print "temp:",temp
            ising1D.system.h = h
            #
            line = []
            #
            for d in range(1,d_max+1):
                line.append( ising1D.canonical.thermal_xxcorrelation(d, temp ))         
            if (h<1.0):
                xi_theory.append( exact(h,temp) )
            else:
                xi_theory.append( exact2(h,temp) )

            data = numpy.log(line[d_min_fit-1:d_max_fit])
            xdata = numpy.arange(d_min_fit,d_max_fit+1)

            p = numpy.polyfit(xdata,data,1)

            axes.plot(xdata,data,"o")
            axes.plot(xdata,p[0]*xdata+p[1],"-")
            
            xi.append( -1.0/p[0] )
    #
    numpy.savetxt( datadir+'/data.dat', numpy.array([ temp_vect,xi ]).transpose() ,fmt='%.18e')

    numpy.savetxt( datadir+'/data2.dat', numpy.array( [temp_vect, xi_theory ] ).transpose() ,fmt='%.18e')
    #
    plot(False)
    #
    return

def plot(append):
    """
    It plots the results of calc reading the file data.dat
    """
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir + '/data.dat')
    data2 = numpy.loadtxt( datadir + '/data2.dat')
    #
    if (append):
        ax = pp.gca()
    else:
        fig = pp.figure()
        ax  = fig.add_subplot(111)
    
    ax.plot(numpy.log10(data[:,0]), numpy.log10(data[:,1]) ,".",ms=8 )
    
    ax.plot(numpy.log10(data2[:,0]), numpy.log10(data2[:,1]) ,"-",lw = 2.0,ms=3 )
    #
    return
