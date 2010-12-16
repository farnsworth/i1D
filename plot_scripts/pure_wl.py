
import sys
sys.path.append('../libs/')
sys.path.append('../')

from all import datadir
import numpy
import math
import ising1D
import utilities
import system
import matplotlib.pyplot as pp
from  time import sleep,time

from matplotlib import rc
rc('text', usetex=True)

def wl( h, L, emin, emax, nbin , nsim , all = True , threshold = 0.2 ):
    #
    #
    ising1D.system.l = L
    ising1D.system.h = h
    ising1D.wl.threshold = threshold
    #
    print"\n"
    print "WL calculation"

    for i in range(nsim):
        if (i==0):
            read=False
            # repeat the calculation if i'm not able to find an initial conf
            j = 0
            while (j<10):
                j = j+1
                random_seed = numpy.random.random_integers(2024, size=8)
                res = ising1D.wl.wl_edos(emin, emax, nbin , random_seed, read)
#                return
                if (res==0): break
        else:
            read=True
            random_seed = numpy.random.random_integers(2024, size=8)
            res = ising1D.wl.wl_edos(emin,emax, nbin , random_seed,read)

 
        
        if (res!=0):
            print "ERRROR: some error in the calculation"
            return
        
        if (i==0):
            xdata = ising1D.wl.logdosf[:,0].copy()
            data  = numpy.zeros( len(xdata) )
            data2 = data.copy()
            delta = ising1D.wl.delta
        
        print "<<< number of bins >>>",len(xdata)
        if (len(ising1D.wl.logdosf[:,0]) != len(xdata) ):
            print "error:",len(xdata),len(ising1D.wl.logdosf[:,0])
            return
        
        res = utilities.HistoFromLog(ising1D.wl.logdosf[:,1],delta)
        data += res
        data2 += res*res
    #
    data = data/numpy.float(nsim)
    data2 = data2/numpy.float(nsim)
    #
    sigma = numpy.sqrt(data2 - data*data)
    sigma = sigma/numpy.sqrt( numpy.float(nsim) )
    #
    tmp = numpy.column_stack( ( xdata , data, sigma ) )
    if (all):
        numpy.savetxt( datadir+'data0.out', tmp, fmt='%.18e' )
        do_plot()
    else:
        numpy.savetxt( datadir+'data1.out', tmp, fmt='%.18e' )
        do2_plot( emin, emax, delta )


    return


def do_plot():
    data = numpy.loadtxt( datadir+'data0.out')
#    data3 = numpy.loadtxt( datadir+'data3.out')

    fig = pp.figure()
    ax = fig.add_subplot(111)

    ax.plot(data[:,0],data[:,1],'-')
    ax.errorbar(data[:,0],data[:,1],yerr = 3.0*data[:,2],fmt=None)

    return

def do2_plot( emin, emax, delta ):
    data_all = select_data_from_all( emin, emax, delta )
    data = numpy.loadtxt( datadir + 'data1.out' )
    
    fig = pp.figure()
    ax = fig.add_subplot(111)

    ax.plot(data[:,0],data[:,1],'.',label="restricted")
    ax.errorbar(data[:,0],data[:,1],yerr = 3.0*data[:,2],fmt=None)

    ax.plot(data_all[:,0],data_all[:,1],'-',label="from all")
    ax.errorbar(data_all[:,0],data_all[:,1],yerr = 3.0*data_all[:,2],fmt=None)
    ax.legend( loc = 0)

    return


def select_data_from_all( emin, emax , deltae ):
    data = numpy.loadtxt( datadir+'data0.out')
    data_new = []
    norm = 0.0
    norm_err = 0.0
    for i in range(len(data[:,0])):
        if ( (data[i,0]<emax) and ( data[i,0]>emin ) ):
            data_new.append( data[i,:] )
            norm = norm + data[i,1]
            norm_err = norm_err + data[i,2]*data[i,2]

    data_new = numpy.array(data_new)
    norm_err = numpy.sqrt(norm_err)
    rel_err = data_new[:,2]/data_new[:,1]
    norm_rel_err = norm_err/norm
    tot_rel_err = numpy.sqrt( norm_rel_err**2 + rel_err**2 )
    data_new[:,1] = data_new[:,1]/(norm*deltae)
    data_new[:,2] = data_new[:,1]*tot_rel_err
    
    return data_new


def wl_calpha(h0,h, l, nsim, enbins, threshold = 0.2):

    ising1D.system.l = l
    ising1D.system.h = h
    ising1D.system.h0 = h0
    ising1D.wl_new.threshold = threshold

    random_seed = numpy.random.random_integers(2024, size=8)
    ising1D.wl_new.wl_ediag(enbins, random_seed )
    delta = ising1D.wl_new.delta
    res = utilities.HistoFromLog(ising1D.wl_new.logdosf[:,1],delta)
    
    data = res.copy()
    data2 = data*data

    for i in range(nsim-1):
        random_seed = numpy.random.random_integers(2024, size=8)
        ising1D.wl_new.wl_ediag(enbins, random_seed )
        res = utilities.HistoFromLog(ising1D.wl_new.logdosf[:,1],delta)
        data = data + res
        data2 = data2 + res*res

    ising1D.quench.en_dist_calc( enbins )
    data_exact = ising1D.quench.dist.copy()


    data = data/numpy.float(nsim)
    data2 = data2/numpy.float(nsim)
    sigma = numpy.sqrt( data2 - data*data )/numpy.sqrt(numpy.float(nsim))

    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.plot(ising1D.wl_new.logdosf[:,0], data , label="wl" )
    ax.plot( data_exact[:,0],data_exact[:,1], "-", label="exact")
    ax.errorbar(ising1D.wl_new.logdosf[:,0],data,yerr = 3.0*sigma,fmt=None)
    ax.legend(loc=0)

    err = numpy.abs(data - data_exact[:,1])/numpy.abs(sigma)

    ax3 = pp.axes( [ 0.5 , 0.2 , 0.3 , 0.3 ] )
    ax3.plot( data_exact[:,0], err, 'go')
#    ax3.set_ylim(0,5)

    return



def wl_calpha_mag (h0 , h, l, nsim, mnbins, threshold = 0.2, accuracy = 1.0e-6):
    #
    ising1D.system.l = l
    ising1D.system.h = h
    ising1D.system.h0 = h0
    ising1D.wl_new.threshold = threshold
    ising1D.wl_new.accuracy = accuracy
    #
    random_seed = numpy.random.random_integers(2024, size=8)
    ising1D.wl_new.wl_mdiag(mnbins, random_seed )
    delta = ising1D.wl_new.delta
#    res = utilities.HistoFromLog(ising1D.wl_new.logdosf[:,1],delta)
    res = ising1D.wl_new.logdosf[:,1]
    
    data = res.copy()
    data2 = data*data

    for i in range(nsim-1):
        random_seed = numpy.random.random_integers(2024, size=8)
        ising1D.wl_new.wl_mdiag(mnbins, random_seed )
#        res = utilities.HistoFromLog(ising1D.wl_new.logdosf[:,1],delta)
        res = ising1D.wl_new.logdosf[:,1]
        data = data + res
        data2 = data2 + res*res

    ising1D.quench.mag_dist_calc( mnbins )
#    data_exact = ising1D.quench.dist.copy()

    data = data/numpy.float(nsim)
    data2 = data2/numpy.float(nsim)
    sigma = numpy.sqrt( data2 - data*data )/numpy.sqrt(numpy.float(nsim))
    
    alpha = -numpy.log( delta*numpy.sum( numpy.exp(data) ) )

    data = data + alpha

    fig = pp.figure()
    ax  = fig.add_subplot(111)
    ax.plot(ising1D.wl_new.logdosf[:,0], data ,".",ms=10)#, label="numerical" )
#    ax.plot( data_exact[:,0],numpy.log(data_exact[:,1]), "-", linewidth=2.0 , label="exact")
#    ax.errorbar(ising1D.wl_new.logdosf[:,0],data,yerr = sigma,fmt=None)
#    leg = ax.legend(loc=0)
#    for t in leg.get_texts():
#        t.set_fontsize(20)

    ax.set_xlabel(r"$M_z$",fontsize=25)
    ax.set_ylabel(r"$\log P_\mathrm{D}$",fontsize=25 )

##    for i in (range(len(sigma))):
##        if (sigma[i] == 0.0 ):
##            print "errore",i

    
#    err = numpy.abs((data - numpy.log(data_exact[:,1]) ))/sigma

#    rel_err = numpy.log( numpy.abs( (data - numpy.log(data_exact[:,1]))/numpy.log(data_exact[:,1]) ) )

##    rel_err =  numpy.abs( data - numpy.log(data_exact[:,1]) )

##    rel_err = numpy.log( numpy.abs( (numpy.exp(data) - data_exact[:,1])/data_exact[:,1]) )

#    ax3 = pp.axes( [ 0.5 , 0.2 , 0.35 , 0.35 ] )
#    ax3.plot( data_exact[:,0], rel_err, 'go')
#    ax3.set_xlim(-0.8,0.8)
#    ax3.set_xticks( [-0.5,0,0.5] )
#    ax3.set_yticks( [-2,-6,-10] )
#    ax3.set_xlabel(r"$M_z$")
#    ax3.set_ylabel(r"$\log \epsilon$")
##    ax3.set_ylim(0,5)

    [i.set_fontsize(15) for i in ax.get_xticklabels() + ax.get_yticklabels()]

##    ax.set_xlim(-1.1)
    ax4 = pp.axes( [ 0.5 , 0.2 , 0.35 , 0.35 ] )
    ax4.plot( ising1D.wl_new.logdosf[1:,0], sigma[1:], 'go')
    ax4.set_xlim(-0.8,0.8)
    ax4.set_xticks( [-0.5,0,0.5 ] )
#    ax4.set_yticks( [0.0,0.01,0.02,0.03] )
##    ax4.set_ylabel(r"$\delta / \sigma$")
    ax4.set_xlabel(r"$M_z$")
    ax4.set_ylabel(r"$\sigma$")

    return
