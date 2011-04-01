
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

def get_last_dat():
    files = os.listdir(datadir);
    dat = []
    [dat.append(int(name[4:6])) for name in files if name.endswith(".dat")]
    if (len(dat)>0):
        return numpy.array(dat).max()
    else:
        return 0

def calc( h0_in, h_in, l ):
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.l = l
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    #h0_in = numpy.array(h0_in)
    #h_in = numpy.array(h_in)
    dmax = 100
    data = []
    #
    fit_fig_can = pp.figure()
    fit_ax_can  = fit_fig_can.add_subplot(111)
    #
    #fit_fig_gge = pp.figure()
    #fit_ax_gge  = fit_fig_gge.add_subplot(111)
    #
    for h0 in h0_in:
        for h in h_in:
            if (numpy.abs(h- h0) < 1.0e-5):
                continue
            ising1D.system.h = h
            ising1D.system.h0 = h0
            xxcorrelation1 = []
            xxcorrelation2 = []
            xxcorrelation3 = []
            #
            energy = system.E0(h0,h,l)
            temperature = ising1D.canonical.find_temperature( energy , 0.000001 )
            #
            #
            #print "effective temperature",temperature
            #
            for d in range(1,dmax):
                xxcorrelation1.append( ising1D.canonical.thermal_xxcorrelation(d,temperature))
                xxcorrelation2.append( ising1D.quench.long_time_xxcorrelation(d) )
            #    xxcorrelation3.append( ising1D.gge.gge_xxcorrelation(d))
            #

            #
            # canonical calculation
            #
            fit_ax_can.plot( range(1,dmax), numpy.log10(numpy.abs(numpy.array(xxcorrelation2))),'.',label = "h0 = %.1f" %h0 )
            pp.draw()
            #
            dmin_fit = int(raw_input("dmin >"))
            dmax_fit = int(raw_input("dmax >"))
            if (dmax_fit > dmax -1):
                print "dmax too big"
                return
            # data to fit
            #
            data1 = numpy.log(numpy.abs(numpy.array(xxcorrelation1[dmin_fit-1:dmax_fit])))
            xdata = numpy.arange(dmin_fit,dmax_fit+1)
            p1 = numpy.polyfit(xdata,data1,1)
            #print p1
            data2 = numpy.log(numpy.abs(numpy.array(xxcorrelation2[
dmin_fit-1:dmax_fit])))
            p2 = numpy.polyfit(xdata,data2,1)
            #print p2
            xdata = numpy.arange(1,dmax)
            fit_ax.plot(xdata,(xdata*p2[0]+p2[1])*numpy.log10(numpy.exp(1)),'-',lw=2.0)
            pp.draw()
            #
            #
            data.append( [h0,h,energy,temperature,numpy.abs(1.0/p1[0]),numpy.abs(1.0/p2[0]) ] )
    #
    fit_ax.legend(loc=0)
    pp.draw()
    #
    index = get_last_dat() + 1
    numpy.savetxt( datadir+'/data'+('%02d' %index)+'.dat', numpy.array(data) ,fmt='%.18e')
    print "Data saved in data"+('%02d' %index)+'.dat'
    #
    plot( index )
    #
    return

def plot( index ):
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir + '/data'+('%02d' %index)+'.dat' )
    print " Data with h = ",data[0,1]
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    #
    ax.plot( data[:,0], numpy.log10(numpy.abs(data[:,4])) ,".",ms=4 ,label="thermal")
    ax.plot( data[:,0], numpy.log10(numpy.abs(data[:,5])) ,".",ms=4 ,label="longtime")
    #
    ax.legend(loc=0)
    #
    fig2 = pp.figure()
    ax2  = fig2.add_subplot(111)
    ax2.plot( numpy.log10(data[:,3]), numpy.log10(numpy.abs(data[:,4])) ,"-",ms=4 ,label="thermal")
    ax2.plot( numpy.log10(data[:,3]), numpy.log10(numpy.abs(data[:,5])) ,".",ms=4 ,label="long time")
    ax2.set_xlim([-2.0,0.0])
    ax2.set_ylim([0.0,5.0])
    #
    #
    #
    return


def append_plot(index):
    #
    print "... Appending the plot ..."
    #
    data = numpy.loadtxt( datadir + '/data'+('%02d' %index)+'.dat' )
    print " Data with h = ",data[0,1]
    #
    if (len(pp.get_fignums())==0):
        print "No figure where append"
        return
    
    ax = pp.gca()

    ax.plot( numpy.log10(data[:,3]), numpy.log10(numpy.abs(data[:,4])) ,"-",ms=4 ,label="thermal")
    ax.plot( numpy.log10(data[:,3]), numpy.log10(numpy.abs(data[:,5])) ,".",ms=4 ,label="long time")

    pp.draw()
