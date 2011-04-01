
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

def calc( h0_in, h_in, l, fstr, dmax ,time=100.0, defaultmin=30):
    """
    h0_in    : is a list of sarting h
    h_in     : is a list of final h
    l        : length of the system
    fstr     : function to use for calculating the xxcorrelation
               examples: "thermal", "gga", "gs","long_time", "time"
    dmax=100 : max value of correlation distance
    """
    print "...Doing calculation..."
    #
    ising1D.system.l = l
    ising1D.system.datadir[:len(datadir)] = datadir
    #
    data = []
    #
    fit_fig = pp.figure()
    fit_ax  = fit_fig.add_subplot(111)
    #
    if (fstr=="thermal"):
        print "thermal calculation"
        function = ising1D.canonical.thermal_xxcorrelation
    elif (fstr=="gge"):
        print "gge calculation"
        function = ising1D.gge.gge_xxcorrelation
    elif (fstr=="gs"):
        print "gs calculation"
        function = ising1D.exact.gs_xxcorr
    elif (fstr=="time"):
        print "time"
        function = ising1D.quench.time_xxcorrelation
    else:
        print "longtime calculation"
        function = ising1D.quench.long_time_xxcorrelation
    #
    for h0 in h0_in:
        for h in h_in:
            #
            ising1D.system.h = h
            ising1D.system.h0 = h0
            xxcorrelation = []
            print "h",h,"h0",h0
            #
            if (numpy.abs(h- h0) < 1.0e-5):
                # becouse the exponential behaviour is present 
                # for very large l
                continue
                energy = system.E0(h0,h,l)
                temperature = 0.0
                if (h<1.0):
                    continue
            else:
                energy = system.E0(h0,h,l)
                temperature = ising1D.canonical.find_temperature( energy , 0.000001 )
            #
            #
            for d in range(1,dmax):
                if (fstr=="thermal"):
                    a = function(d,temperature)
                    xxcorrelation.append( a )
                elif ( fstr =="time"):
                    xxcorrelation.append( function(time,d) )
                else:
                    xxcorrelation.append( function(d))
            #
            temp = numpy.log10(numpy.abs(numpy.array(xxcorrelation)))
            fit_ax.plot( range(1,dmax), temp ,'-',label = r"$h_0 = %.1f$" %h0,lw=2.0 )
            
            fit_ax.set_ylim(top=0,bottom=temp.min())
            pp.draw()
            #
            try:
                dmin_fit = int(raw_input("dmin >"))
            except:
                dmin_fit = defaultmin

            try:
                dmax_fit = int(raw_input("dmax >"))
            except:
                dmax_fit = dmax - 1

            if (dmax_fit > dmax -1):
                dmax_fit = dmax - 1
    
            data_corr = numpy.log(numpy.abs(numpy.array(xxcorrelation[dmin_fit-1:dmax_fit])))
            xdata = numpy.arange(dmin_fit,dmax_fit+1)
            p = numpy.polyfit(xdata,data_corr,1)

            xdata = numpy.arange(1,dmax)
            tmp = -1.0/p[0]
            fit_ax.plot(xdata,numpy.log10(numpy.exp(1.0))*(xdata*p[0]+p[1]),'--',lw=1.0,label=r"$\xi= %.1f$" %tmp)
            pp.draw()
            #
            #
            data.append( [h0,h,energy,temperature,-1.0/p[0]] )
            print "calculated value",-(1.0/p[0])
    #
    fit_ax.legend(loc=0)
    fit_ax.set_xlabel(r'$l$')
    fit_ax.set_ylabel(r'$\log ( \rho^{xx}_l ) $')
    pp.legend(loc=0)
    pp.draw()
    #
    index = get_last_dat(fstr)+1
    numpy.savetxt( datadir+'/data_'+fstr+'%02d.dat' %index, numpy.array(data) ,fmt='%.18e')
    print "Data saved in data_"+fstr+'%02d.dat' %index
    #
    if len(h0_in)*len(h_in)>1:
        plot( [fstr+'%02d' %index] )
    #
    if (h_in[0]>1.0):
        print "zero temperature xi",1.0/numpy.log(h_in[0])
    #
    return index

def plot( functions ):
    #
    print "... Doing the plot ..."
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    #
    for fstr in functions:
        data = numpy.loadtxt( datadir + '/data_'+fstr+'.dat' )
        print "final h =",data[0,1]
        #
        #if ( fstr.find("thermal")>=0):
        #    print "funzia"
        #    ax.plot( 0.5*data[:,3],0.5*numpy.abs(data[:,4]) ,".",ms=4 ,label=fstr)
        #else:
        ax.plot( data[:,3], numpy.abs(data[:,4]) ,".",ms=4 ,label=fstr)
    #
    ax.legend(loc=0)
    #
    #
    return



######## plot 1 ########
# generate plot fig5 Rossini et al.
def generate_plot(function):
    calc( [1.1,1.2,1.5,2.0], [1.25], 1000, function,100)

    index = calc( numpy.arange(1.05,1.8,0.05), [1.25], 1000, function,dmax=150)
    raw_input("Select the figure and then Press Enter ...")
    fig = pp.gcf()
    plotsmall(fig,function+'%02d' %index,1.25)

    return

def plotsmall(fig,fstr,h):
    xigs = 1.0/numpy.log(h)
    ax2 = fig.add_axes( (0.2,0.18,0.35,0.3) )
    ax2.set_xlabel(r"$h_0$")
    ax2.set_ylabel(r"$\xi$")
    data = numpy.loadtxt( datadir + '/data_'+fstr+'.dat' )
    ax2.plot( data[:,0], numpy.abs(data[:,4]) ,".",ms=4 ,label=fstr)
    ax2.plot([h],[xigs],"*",ms=6 )
    pp.draw()



######## plot 2 ########
# generate plot fig7 Rossini at al.
def generate_plot2():

    fig = pp.figure()
    ax = fig.add_subplot(111)
    ax.set_title("prova")
    
    function="time"
    index = calc( numpy.arange(1.05,1.25,0.05), [1.25], 1000, function,dmax=150)
    plotxi(ax,function+'%02d' %index)
    pp.draw()

    index = calc( numpy.arange(1.25,1.8,0.05), [1.25], 1000, function,dmax=150)
    raw_input("Select the figure and then Press Enter ...")
    ax = pp.gca()
    plotxi(ax,function+'%02d' %index)
    pp.draw()

    function="gge"
    index = calc( numpy.arange(1.05,1.25,0.05), [1.25], 1000, function,dmax=150)
    raw_input("Select the figure and then Press Enter ...")
    ax = pp.gca()
    plotxi(ax,function+'%02d' %index)
    pp.draw()

    index = calc( numpy.arange(1.25,1.8,0.05), [1.25], 1000, function,dmax=150)
    raw_input("Select the figure and then Press Enter ...")
    ax = pp.gca()
    plotxi(ax,function+'%02d' %index)
    pp.draw()

    function="thermal"
    index = calc( numpy.arange(1.25,1.8,0.05), [1.25], 1000, function,dmax=150)
    raw_input("Select the figure and then Press Enter ...")
    ax = pp.gca()
    plotxi(ax,function+'%02d' %index)

    ax.set_xlabel(r"$T_{\mathrm{eff}}$")
    ax.set_ylabel(r"$\xi$")
    ax.legend(loc=0)
    pp.draw()

    return

def plotxi(ax,fstr):
    data = numpy.loadtxt( datadir + '/data_'+fstr+'.dat' )
    ax.plot( data[:,3], numpy.abs(data[:,4]) ,".",ms=4 ,label=fstr)



# trial with revivals
def generate_plot_temp():
    calc( [1.1,1.2,1.5,2.0], [1.25], 1000, "time",250,time=300.0,defaultmin=120)

    index = calc( numpy.arange(1.05,1.8,0.05), [1.25], 1000, "time",dmax=250,time=300.0,defaultmin=120)
    #raw_input("Select the figure and then Press Enter ...")
    #fig = pp.gcf()
    plotsmall2('time'+'%02d' %index)

    index = calc( numpy.arange(1.25,1.8,0.02), [1.25], 1000, "thermal",dmax=150,defaultmin=80)
    
    raw_input("Select the figure and then Press Enter ...")
    ax = pp.gca()
    
    h = 1.25
    xigs = 1.0/numpy.log(h)

    data = numpy.loadtxt( datadir + '/data_thermal'+'%02d' %index+'.dat' )
    ax.plot( data[:,3], (data[:,4]*xigs)/(xigs-data[:,4]) ,".",ms=4 ,label='thermal')
    pp.draw()
    
    return

def plotsmall2(fstr):
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    #ax2 = fig.add_axes( (0.2,0.18,0.35,0.3) )
    ax.set_xlabel(r"$h_0$")
    ax.set_ylabel(r"$\xi$")
    data = numpy.loadtxt( datadir + '/data_'+fstr+'.dat' )
    ax.plot( data[:,0], data[:,4] ,".",ms=4 ,label=fstr)
    #ax2.plot([h],[xigs],"*",ms=6 )


    fig2 = pp.figure()
    ax2  = fig2.add_subplot(111)
    ax2.plot( data[:,3], numpy.abs(data[:,4]) ,".",ms=4 ,label=fstr)
    pp.draw()

def plotsmall3():
    h = 1.25
    xigs = 1.0/numpy.log(h)

    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.set_xlabel(r"$T_{eff}$")
    ax.set_ylabel(r"$\xi$")
    data = numpy.loadtxt( datadir + '/data_time26.dat' )
    data2 = numpy.loadtxt( datadir + '/data_thermal04.dat' )
    ax.plot( data[:,3], data[:,4] ,".",ms=4 ,label='time')
    ax.plot( data2[:,3], (data2[:,4]*xigs)/(xigs-data2[:,4]) ,".",ms=4 ,label='thermal')
    ax.legend(loc=0)
    pp.draw()
