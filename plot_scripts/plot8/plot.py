
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

def calc( h0_in, h_in , fit_from_value, l ):
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
    xdata = numpy.arange(fit_from_value+1,dmax)
    #
    fit_fig_can = pp.figure()
    fit_ax_can  = fit_fig_can.add_subplot(111)
    fit_ax_can.set_title("canonical")
    #
    fit_fig_time = pp.figure()
    fit_ax_time  = fit_fig_time.add_subplot(111)
    fit_ax_time.set_title("longtime")
    #
    for h0 in h0_in:
        for h in h_in:
            if (numpy.abs(h- h0) < 1.0e-5):
                continue
            ising1D.system.h = h
            ising1D.system.h0 = h0
            xxcorrelation1 = []
            xxcorrelation2 = []
            #
            energy = system.E0(h0,h,l)
            temperature = ising1D.canonical.find_temperature( energy , 0.000001 )
            #
            #
            #print "effective temperature",temperature
            #
            for d in range(1,dmax):
                xxcorrelation1.append( ising1D.canonical.thermal_xxcorrelation(d, temperature))
                xxcorrelation2.append(ising1D.quench.long_time_xxcorrelation( d ) )
            # data to fit
            #
            data1 = numpy.log(numpy.abs(numpy.array(xxcorrelation1[fit_from_value:])))
            p1 = numpy.polyfit(xdata,data1,1)
            #print p1
            data2 = numpy.log(numpy.abs(numpy.array(xxcorrelation2[fit_from_value:])))
            p2 = numpy.polyfit(xdata,data2,1)
            #print p2
            fit_ax_can.plot(xdata,(xdata*p1[0]+p1[1])*numpy.log10(numpy.exp(1)),'-')
            fit_ax_can.plot( range(1,dmax), numpy.log10(numpy.abs(numpy.array(xxcorrelation1))),'.',label = "h0 = %.2f" %h0 )
            #
            fit_ax_time.plot(xdata,(xdata*p2[0]+p2[1])*numpy.log10(numpy.exp(1)),'-')
            fit_ax_time.plot( range(1,dmax), numpy.log10(numpy.abs(numpy.array(xxcorrelation2))),'.',label = "h0 = %.2f" %h0 )
            #
            #
            data.append( [h0,h,energy,temperature,numpy.abs(1.0/p1[0]),numpy.abs(1.0/p2[0]) ] )
    #
    fit_ax_can.legend(loc=0)
    fit_ax_time.legend(loc=0)
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
    #ax.plot( data[:,0], numpy.log10(numpy.abs(data[:,5])) ,".",ms=4 ,label="longtime")
    #
    ax.legend(loc=0)
    #
    fig2 = pp.figure()
    ax2  = fig2.add_subplot(111)
    ax2.plot( numpy.log10(data[:,3]), numpy.log10(numpy.abs(data[:,4])) ,"-",ms=4 ,label="thermal")
    #ax2.plot( numpy.log10(data[:,3]), numpy.log10(numpy.abs(data[:,5])) ,".",ms=4 ,label="long time")
    ax2.set_xlim([-2.0,0.0])
    ax2.set_ylim([0.0,5.0])
    ax2.legend(loc=0)
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
