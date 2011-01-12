
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

def calc( h0, h, l, deltae , nsim, mnbins, threshold = 0.2, accuracy = 1.0e-6):
    #
    print "...Doing calculation..."
    #
    # ... setting parameters
    ising1D.system.l = l
    ising1D.system.h0 = h0
    ising1D.system.h = h
    ising1D.system.datadir[:len(datadir)] = datadir
    ising1D.wl_new.threshold = threshold
    ising1D.wl_new.accuracy = accuracy
    #
    meanE = system.E0(h0,h,l)
    #
    # ... first calculation
    #
    # ... random seed
    random_seed = numpy.random.random_integers(2024, size=8)
    out = ising1D.wl.wl_mdos(meanE,deltae,mnbins, random_seed,False)
    if (out != 0):
      print "some error in the calcualtion!"
  
#    delta = ising1D.wl_new.delta
#    res = utilities.HistoFromLog(ising1D.wl_new.logdosf[:,1],delta)
    res = ising1D.wl.logdosf[:,1]
    #
    xdata = ising1D.wl.logdosf[:,0].copy()
    data = res.copy()
    data2 = data*data
    #
    for i in range(nsim-1):
        random_seed = numpy.random.random_integers(2024, size=8)
        out = ising1D.wl.wl_mdos( meanE,deltae,mnbins, random_seed,False)
#        res = utilities.HistoFromLog(ising1D.wl_new.logdosf[:,1],delta)
        res = ising1D.wl.logdosf[:,1].copy()
        if (len(res) != len(xdata)):
            print "Error: different length of results"
            return
        data = data + res
        data2 = data2 + res*res
    #
    # ... compute average and sigma
    data = data/numpy.float(nsim)
    data2 = data2/numpy.float(nsim)
    sigma = numpy.sqrt( data2 - data*data )/numpy.sqrt(numpy.float(nsim))
    # ... exact calculation
    ising1D.quench.sigma_mag_quench2(meanE, deltae, mnbins )
    data_exact = []
    for i in range(len(ising1D.quench.dist[:,0])):
        if (ising1D.quench.dist[i,1] > 0.0):
            data_exact.append( ising1D.quench.dist[i,:] )
    #
    # save data for plot
    numpy.savetxt( datadir+'/data_num.dat',numpy.column_stack( (ising1D.wl.logdosf[:,0],data,sigma) ) , fmt='%.18e')
    #
    numpy.savetxt( datadir+'/data_exact.dat', data_exact,fmt='%.18e')
    #
    if (nsim>1):
        plot(True)
    else:
        plot(False)
    #
    return

def plot(sigma = False):
    #
    print "... Doing the plot ..."
    #
    data = numpy.loadtxt( datadir+'/data_num.dat')
    data_exact = numpy.loadtxt( datadir+'/data_exact.dat')
    #
    # ... normalization
    deltas = []
    deltas = [ data[i+1,0]-data[i,0] for i in range(len(data[:,0])-1) ]
    delta = numpy.min(deltas)
    #
    alpha = -numpy.log( delta*numpy.sum( numpy.exp(data[:,1]) ) )
    # ... I can't compute the error over alpha because exp(data^2) is too big
    #sigma_alpha = numpy.sqrt( (data[:,2]*data[:,2])*numpy.exp(data[:,1]*data[:,1]))
    #sigma_alpha = sigma_alpha/numpy.sum( numpy.exp(data[:,1]) )
    #
    data[:,1] = data[:,1] + alpha
    #data[:,2] = numpy.sqrt(data[:,2]*data[:,2] + sigma_alpha*sigma_alpha)
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.plot(data[:,0], data[:,1] ,".",ms=10, label="numerical" )
    ax.plot( data_exact[:,0],numpy.log(data_exact[:,1]), "-", linewidth=2.0 , label="exact")
    if (sigma):
        ax.errorbar( data[:,0],data[:,1],yerr = data[:,2] )#,fmt=None)
    leg = ax.legend(loc=1)
    #
    [t.set_fontsize(15) for t in leg.get_texts()]
    ax.set_xlabel(r"$M_z$",fontsize=25)
    ax.set_ylabel(r"$\log P_\mathrm{mc}$",fontsize=25 )
    [i.set_fontsize(15) for i in ax.get_xticklabels() + ax.get_yticklabels()]
##
##    for i in (range(len(sigma))):
##        if (sigma[i] == 0.0 ):
##            print "errore",i
##
    # skip the first data because has zero sigma
    if (not sigma):
        err = numpy.abs( ( data[1:,1] - numpy.log(data_exact[1:,1]) )/data[1:,1] )
        rel_err = numpy.log( err )
    else:
        rel_err = numpy.abs( ( data[1:,1] - numpy.log(data_exact[1:,1]) )/data[1:,2] )


    ax3 = pp.axes( [ 0.35 , 0.2 , 0.35 , 0.35 ] )
    ax3.plot( data_exact[1:,0], rel_err, 'go')
    ax3.set_xlabel(r"$M_z$")
    if (not sigma):
        ax3.set_ylabel(r"$\log \epsilon$")
    else:
        ax3.set_ylabel(r"$\epsilon/\sigma$")

#    ax3.set_xlim(-0.8,0.8)
#    ax3.set_xticks( [-0.5,0,0.5] )
#    ax3.set_yticks( [-2,-6,-10] )
##    ax3.set_ylim(0,5)
##    ax.set_xlim(-1.1)

# ... absolute error
#    ax4 = pp.axes( [ 0.5 , 0.2 , 0.35 , 0.35 ] )
#    ax4.plot( ising1D.wl_new.logdosf[1:,0], sigma[1:], 'go')
#    ax4.set_xlim(-0.8,0.8)
#    ax4.set_xticks( [-0.5,0,0.5 ] )
#    ax4.set_yticks( [0.0,0.01,0.02,0.03] )
#    ax4.set_xlabel(r"$M_z$")
#    ax4.set_ylabel(r"$\sigma$")

    return
