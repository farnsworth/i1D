
import sys
import os
sys.path.append('../../libs/')
sys.path.append('../../')

datadir=os.getcwd()

import numpy
import math
import ising1D
import utilities
import system
import matplotlib.pyplot as pp
from  time import sleep,time

from matplotlib import rc
rc('text', usetex=True)


def calc(h0 , h, l, nsim, mnbins, threshold = 0.2, accuracy = 1.0e-6):
    #
    # ... setting parameters
    ising1D.system.l = l
    ising1D.system.h = h
    ising1D.system.h0 = h0
    ising1D.system.datadir = datadir
    ising1D.wl_new.threshold = threshold
    ising1D.wl_new.accuracy = accuracy
    #
    # ... random seed
    random_seed = numpy.random.random_integers(2024, size=8)
    #
    # ... first calculation
    ising1D.wl_new.wl_mdiag(mnbins, random_seed )
    delta = ising1D.wl_new.delta
#    res = utilities.HistoFromLog(ising1D.wl_new.logdosf[:,1],delta)
    res = ising1D.wl_new.logdosf[:,1]
    #
    data = res.copy()
    data2 = data*data

    for i in range(nsim-1):
        random_seed = numpy.random.random_integers(2024, size=8)
        ising1D.wl_new.wl_mdiag(mnbins, random_seed )
#        res = utilities.HistoFromLog(ising1D.wl_new.logdosf[:,1],delta)
        res = ising1D.wl_new.logdosf[:,1]
        data = data + res
        data2 = data2 + res*res
    #
    # ... compute average and sigma
    data = data/numpy.float(nsim)
    data2 = data2/numpy.float(nsim)
    sigma = numpy.sqrt( data2 - data*data )/numpy.sqrt(numpy.float(nsim))
    #
    # ... normalization
    alpha = -numpy.log( delta*numpy.sum( numpy.exp(data) ) )
    data = data + alpha
    #
    # ... exact calculation
    ising1D.quench.mag_dist_calc( mnbins )
    data_exact = ising1D.quench.dist.copy()
    #
    # save data for plot
    numpy.savetxt( datadir+'data_num.dat',numpy.column_stack( (ising1D.wl_new.logdosf[:,0],data,sigma) ) , fmt='%.18e')
    #
    numpy.savetxt( datadir+'data_exact.dat', data_exact)
    #
    return


def plot():
    #
    data = numpy.loadtxt( datadir+'data_num.dat')
    data_exact = numpy.loadtxt( datadir+'data_exact.dat')
    #
    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.plot(data[:,0], data[:,1] ,".",ms=10, label="numerical" )
    ax.plot( data_exact[:,0],numpy.log(data_exact[:,1]), "-", linewidth=2.0 , label="exact")
    ax.errorbar(ising1D.wl_new.logdosf[:,0],data,yerr = sigma,fmt=None)
    leg = ax.legend(loc=0)
    #
    [t.set_fontsize(20) for t in leg.get_texts()]
    ax.set_xlabel(r"$M_z$",fontsize=25)
    ax.set_ylabel(r"$\log P_\mathrm{D}$",fontsize=25 )
    [i.set_fontsize(15) for i in ax.get_xticklabels() + ax.get_yticklabels()]
##
##    for i in (range(len(sigma))):
##        if (sigma[i] == 0.0 ):
##            print "errore",i
##  
    err = numpy.abs((data[:,1] - numpy.log(data_exact[:,1]) ))/data[:,2]
    rel_err = numpy.log( err )

    ax3 = pp.axes( [ 0.5 , 0.2 , 0.35 , 0.35 ] )
    ax3.plot( data_exact[:,0], rel_err, 'go')
    ax3.set_xlabel(r"$M_z$")
    ax3.set_ylabel(r"$\log \epsilon$")

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
