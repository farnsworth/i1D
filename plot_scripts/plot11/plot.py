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
import scipy.special
#from  time import sleep,time

from matplotlib import rc
rc('text', usetex=True)

def exactPfeuty(h,d):
    lamb = 1/h
    temp = ((1-lamb*lamb)**(-0.25))*((d*numpy.pi)**(-0.5))*(lamb**d)
    return temp
    
def exactSachdev(h,temp,d):
    c = 2.0*numpy.sqrt(h)
    m = (1.0-h)/(2.0*h)
    delta = 2.0*numpy.abs(1-h)
    xi = numpy.exp(delta/temp)*numpy.sqrt(numpy.pi*c/(2.0*delta*temp))
    Rt = numpy.exp(-d/xi)
    x = delta*d/c
    K0 = scipy.special.k0(x)
    #print Rt,K0

    return Rt*K0

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
    xxcorrelations = []
    xxcorrelations2 = []
    #
    xxcorrelations.append(range(1,d_max+1))
    xxcorrelations2.append(range(1,d_max+1))
    #
    #
    d_max_fit = d_max
    d_min_fit = 50
    #
    #
    for h in h_vect:
        for temp in temp_vect:

            h_limit = 0.0

            if (h < 1.0):
                h_limit = (1.0-h*h)**(0.25)
            #
            iszerotemp = temp < 1.0e-6
            ising1D.system.h = h
            ising1D.system.h0 = h
            line = []
            line2 = []
            for d in range(1,d_max+1):
                line.append( ising1D.canonical.thermal_xxcorrelation(d, temp ))
                if (iszerotemp):
                    line2.append( exactPfeuty(h,d) )
                else:
                    line2.append( exactSachdev(h,temp,d) )

            xxcorrelations.append(line)
#            if (iszerotemp):
            xxcorrelations2.append(line2)

            data = numpy.log(line[d_min_fit-1:d_max_fit])
            xdata = numpy.arange(d_min_fit,d_max_fit+1)

#        print len(data)
#        print len(xdata)

            p = numpy.polyfit(xdata,data,1)

            if (iszerotemp):
                print "h",h,"length, theoretical",-numpy.log(numpy.log(h)),"from fit",-numpy.log(-p[0])

    #
    numpy.savetxt( datadir+'/data.dat', numpy.array(xxcorrelations).transpose() ,fmt='%.18e')
    if (len(xxcorrelations2)>0):
        numpy.savetxt( datadir+'/data2.dat', numpy.array(xxcorrelations2).transpose() ,fmt='%.18e')
    #
    plot(len(xxcorrelations2)>0)
    #
    return

def plot(alsoexact):
    """
    It plots the results of calc reading the file data.dat
    """
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir + '/data.dat')
    if (alsoexact):
        data2 = numpy.loadtxt( datadir + '/data2.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    for i in range(len(data[1,:])-1):
        ax.plot(data[:,0], numpy.log(data[:,i+1]) ,".",ms=8 )
        
    if (alsoexact):
        for i in range(len(data2[1,:])-1):
            ax.plot(data2[:,0], numpy.log(data2[:,i+1]) ,"-",lw = 2.0,ms=3 )
    #
    return
