
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

def get_last_dat(root):
    files = os.listdir(datadir);
    dat = []
    [dat.append(int(name[len(name)-6:len(name)-4])) for name in files if (name.endswith(".dat") and name.startswith("data_"+root) )]
    if (len(dat)>0):
        return numpy.array(dat).max()
    else:
        return 0

def calc( h0_in, h_in,dmax=100, l=1000 ):
    """
    h0_in    : is a list of sarting h
    h_in     : is a list of final h
    l = 1000       : length of the system
    """
    print "...Doing calculation..."
    #
    ising1D.system.l = l
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    datalt = [ range(1,dmax) ]
    datagge = [ range(1,dmax)]
    #
    for h0 in h0_in:
        for h in h_in:
            #
            ising1D.system.h = h
            ising1D.system.h0 = h0
            xxcorrelation1 = []
            xxcorrelation2 = []
            #
            for d in range(1,dmax):
                xxcorrelation1.append( ising1D.quench.long_time_xxcorrelation(d) )
                xxcorrelation2.append( ising1D.gge.gge_xxcorrelation(d) )
            #
            datalt.append(xxcorrelation1)
            datagge.append(xxcorrelation2)
            #
    #
    #index = get_last_dat(fstr)+1
    numpy.savetxt( datadir+'/datalt.dat', numpy.array(datalt).transpose() ,fmt='%.18e')
    numpy.savetxt( datadir+'/datagge.dat', numpy.array(datagge).transpose() ,fmt='%.18e')
    #
    plot()
    #
    return

def plot( ):
    #
    print "... Doing the plot ..."
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)

    datalt = numpy.loadtxt( datadir + '/datalt.dat' )
    datagge = numpy.loadtxt( datadir + '/datagge.dat' )

    for i in range(1,len(datalt[1,:])):
        ax.plot( datalt[:,0], numpy.abs(datalt[:,i]-datagge[:,i]) ,".",ms=4 ,label="long time")
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    #
    dmin_fit = 50
    xdata = datalt[dmin_fit:,0]
    #
    for i in range(1,len(datalt[1,:])):
        data = numpy.log(numpy.abs(datalt[dmin_fit:,i]/datagge[dmin_fit:,i]))
        data2 = numpy.log(numpy.abs(datagge[dmin_fit:,i]))
        p = numpy.polyfit(xdata,data,1)
        p2 = numpy.polyfit(xdata,data2,1)
        
        ax.plot( datalt[:,0], numpy.log(numpy.abs(datalt[:,i]/datagge[:,i])) ,".",ms=4 ,label="$\epsilon = %.2e$" %(p[0]*p2[0]))
        ax.plot( datalt[:,0], datalt[:,0]*p[0]+p[1])
    #
    pp.legend(loc=0)
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    #
    for i in range(1,len(datalt[1,:])):
        ax.plot( datalt[:,0], numpy.abs((datalt[:,i]-datagge[:,i])/datalt[:,i]) ,".",ms=4 ,label="long time")
    #
    #
    return
