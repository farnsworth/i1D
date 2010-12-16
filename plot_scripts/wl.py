
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

def Conf( h0, h, L, deltae, nbin, nsim, threshold, accuracy = 1.0e-6 ):
    #
    #
    ising1D.system.l = L
    ising1D.system.h = h
    ising1D.system.h0 = h0
    ising1D.wl.threshold = threshold
    ising1D.wl.accuracy = accuracy
    #
    meanE = system.E0(h0,h,L)
    #
    print "Exact calculation"
    ising1D.quench.sigma_mag_quench2( meanE, deltae, nbin)
    #
    delta = ising1D.quench.dist[1,0]-ising1D.quench.dist[0,0]
    #
    sleep(5)
    print"\n"
    print "WL calculation"

    for i in range(nsim):
        j = 0
        if (i==0):
            read=False
            # repeat the calculation if i'm not able to find an initial conf
            while (j<10):
                random_seed = numpy.random.random_integers(2024, size=8)
                res = ising1D.wl.wl_mdos(meanE,deltae,nbin,random_seed,read)
                if (res==0): break
        else:
            read=True
            random_seed = numpy.random.random_integers(2024, size=8)
            res = ising1D.wl.wl_mdos(meanE,deltae,nbin,random_seed,read)

 
        
        if (res!=0):
            print "ERRROR: some error in the calculation"
            return
        
        if (i==0):
            xdata = ising1D.wl.logdosf[:,0].copy()
            data  = numpy.zeros( len(xdata) )
            data2 = data.copy()
        
        if (len(ising1D.wl.logdosf[:,0]) != len(xdata) ):
            
            print "error:",len(xdata),len(ising1D.wl.logdosf[:,0])
            return

#        res = utilities.HistoFromLog(ising1D.wl.logdosf[:,1],delta)
        res = ising1D.wl.logdosf[:,1]
        data += res
        data2 += res*res
    #
    data = data/numpy.float(nsim)
    data2 = data2/numpy.float(nsim)
    #
    alpha = - numpy.log( delta*numpy.sum( numpy.exp(data) ) )
    #
    sigma = numpy.sqrt(data2 - data*data)
    sigma = sigma/numpy.sqrt( numpy.float(nsim) )
    #
    tmp = numpy.column_stack( ( xdata , data+alpha, sigma ) )
    #
    tmp2 = []
    for i in range(len(ising1D.quench.dist[:,1])):
        if (ising1D.quench.dist[i,1]>0.0):
            tmp2.append( [ ising1D.quench.dist[i,0] , numpy.log( ising1D.quench.dist[i,1] ) ] )


#    tmp2 = numpy.column_stack( ( ising1D.quench.dist[:,0] , numpy.log( ising1D.quench.dist[:,1]) )
    numpy.savetxt( datadir+'data1.out', tmp2, fmt='%.18e' )
    numpy.savetxt( datadir+'data2.out', tmp, fmt='%.18e' )

    Conf_plot()

    return


def Conf_plot():
    data1 = numpy.loadtxt( datadir+'data1.out')
    data2 = numpy.loadtxt( datadir+'data2.out')
#    data3 = numpy.loadtxt( datadir+'data3.out')

    fig = pp.figure()
    ax = fig.add_subplot(111)

    ax.plot(data1[:,0],data1[:,1],'-',label="exact",linewidth=2.0)
    ax.plot(data2[:,0],data2[:,1],'.',ms=10,label="Wang-Landau")
#    ax.plot(data3[:,0],data3[:,1],'.')

#    ax.errorbar(data2[:,0],data2[:,1],yerr = data2[:,2],fmt=None,ecolor="k")
    ax.set_xlabel(r"$M_z$",fontsize=20)
    ax.set_ylabel(r"$\log ( P_\mathrm{mc} )$",fontsize=20 )
    ax.set_xticks([0,0.2,0.4,0.6])
    [i.set_fontsize(18) for i in ax.get_xticklabels() + ax.get_yticklabels()]

    delta = data1[1,0]-data1[0,0]

    rel_err = []
    for i in range(len(data2[:,0])):
        for j in range(len(data1[:,0])):
            if (numpy.abs(data2[i,0] - data1[j,0])<delta/2):
                rel_err.append( ( data1[j,1] - data2[i,1])/data2[i,2] )

    ax3 = pp.axes([0.3, 0.2, 0.35, 0.35])
    ax3.set_xlabel(r"$M_z$",fontsize=15)
    ax3.set_ylabel(r"$\frac{abs( \log P_{\mathrm{mc}} - \log P_{\mathrm{mc}}^{\mathrm{wl}})}{ \sigma}$",fontsize=15)
    ax3.plot( data2[:,0],rel_err,'go')
    ax3.set_yticks([-0.5,0.0,0.5])
    ax3.set_xticks([0,0.2,0.4,0.6])

    leg = ax.legend( loc=1)
    for t in leg.get_texts():
        t.set_fontsize(15)

    return


def Confronto_gaussiana(h0, h, L, deltae, nbin):
    #
    #
    ising1D.system.l = L
    ising1D.system.h = h
    ising1D.system.h0 = h0
    #
    meanE = system.E0(h0,h,L)
    #
    print "Exact calculation"
    ising1D.quench.sigma_mag_quench2( meanE, deltae, nbin)
    #
    print "Gaussian calculation 1"
    ising1D.wl_new.exact( meanE, deltae, nbin)
    data1 = ising1D.wl_new.dist.copy()
    #
    print "Gaussian calculation 2"
    ising1D.wl_new.exact( meanE, deltae/2.0, nbin)
    data2 = ising1D.wl_new.dist.copy()
    #
    #
    tmp = numpy.column_stack(( ising1D.quench.dist , data1, data2 ))
    numpy.savetxt( datadir+'data.out', tmp, fmt='%.18e' )

    Gauss_plot()


    return

def Gauss_plot():
    #
    data = numpy.loadtxt( datadir+'data.out')
    #
    fig = pp.figure()
    ax = fig.add_subplot(111)
    #
    ax.plot(data[:,0],data[:,1],'b-')
    ax.plot(data[:,2],data[:,3],'g-')
    ax.plot(data[:,4],data[:,5],'r-')

    return


def under():
    delta = 1.0
    npoint = 1000
    data = 3.0*(2*numpy.arange(npoint) - npoint)/numpy.float(npoint)
    data1 = []
    data2 = gaus(data,0.0,1.0)
    data3 = gaus(data,0.0,0.5)

    for i in data:
        if (abs(i)<delta):
            data1.append(0.5)
        else:
            data1.append(0.0)

    fig = pp.figure()
    ax = fig.add_subplot(111)
    ax.plot(data,data1)
    ax.plot(data,data2)
    ax.plot(data,data3)
    return
    
        
def gaus(x,mu,sigma):
    gaus = (x-mu)*(x-mu)/(2*sigma*sigma)
    gaus = 1.0/(sigma * numpy.sqrt(2.0*numpy.pi) ) * numpy.exp(-gaus)
    return gaus


##################
# microcanonical #
##################


def Conf2( h0, h, L, deltae, nbin, nsim, threshold ):
    #
    #
    ising1D.system.l = L
    ising1D.system.h = h
    ising1D.system.h0 = h0
    ising1D.wl.threshold = threshold
    ising1D.wl_new.threshold = threshold
    #
    meanE = system.E0(h0,h,L)
    #
    print "Exact calculation"
    #ising1D.wl_new.exact( meanE, deltae, nbin)
    #
    ising1D.quench.sigma_mag_quench2( meanE, deltae, nbin )
    #
#    ising1D.quench.sigma_mag_quench2( meanE, 2.0*deltae, nbin)
    #
    delta = ising1D.wl_new.dist[1,0]-ising1D.wl_new.dist[0,0]
    #
    print"\n"
    print "WL calculation"

    for i in range(nsim):
        j = 0
        if (i==0):
            read=False
            # repeat the calculation if i'm not able to find an initial conf
            while (j<10):
                j = j+1
                random_seed = numpy.random.random_integers(2024, size=8)
                res = ising1D.wl_new.wl_mdos(meanE,deltae,nbin,random_seed,read)
                if (res==0): break
                if (res==-1):
                    print "Warning: ",res
                    return
        else:
            read=True
            random_seed = numpy.random.random_integers(2024, size=8)
            res = ising1D.wl_new.wl_mdos(meanE,deltae,nbin,random_seed,read)

 
        
        if (res!=0):
            print "ERROR: some error in the calculation (python message)"
            return
        
        if (i==0):
            xdata = ising1D.wl_new.logdosf[:,0].copy()
            data  = numpy.zeros( len(xdata) )
            data2 = data.copy()
        
        if (len(ising1D.wl_new.logdosf[:,0]) != len(xdata) ):
            
            print "error:",len(xdata),len(ising1D.wl_new.logdosf[:,0])
            return

#        res = utilities.HistoFromLog(ising1D.wl_new.logdosf[:,1],delta)
        res = ising1D.wl_new.logdosf[:,1]
        data += res
        data2 += res*res
    #
    data = data/numpy.float(nsim)
    data2 = data2/numpy.float(nsim)
    #
    sigma = numpy.sqrt(data2 - data*data)
    sigma = sigma/numpy.sqrt( numpy.float(nsim) )
    #
    alpha = -numpy.log( numpy.sum( numpy.exp(data) )    )
    #
    tmp = numpy.column_stack( ( xdata , data+alpha, sigma ) )
#    numpy.savetxt( datadir+'data1.out', ising1D.wl_new.dist, fmt='%.18e' )
    numpy.savetxt( datadir+'data1.out', ising1D.quench.dist , fmt='%.18e' )
    numpy.savetxt( datadir+'data2.out', tmp, fmt='%.18e' )
#    numpy.savetxt( datadir+'data3.out', ising1D.quench.dist, fmt='%.18e' )

    Conf_plot()

    return

def Conf_numerical( h0, h, L, deltae, nbin, nsim, threshold, num_function ):
    # ... stampare tempo impiegato ...
    #
    if (num_function == 0):
        function = ising1D.wl
    else:
        function = ising1D.wl_new

    ising1D.system.l = L
    ising1D.system.h = h
    ising1D.system.h0 = h0
    function.threshold = threshold
    #
    meanE = system.E0(h0,h,L)
    #
    print"\n"
    print "WL calculation"

    t1 = time()

    for i in range(nsim):
        j = 0
        if (i==0):
            read=False
            # repeat the calculation if i'm not able to find an initial conf
            while (j<10):
                j = j+1
                random_seed = numpy.random.random_integers(2024, size=8)
                res = function.wl_mdos(meanE,deltae,nbin,random_seed,read)
                if (res==0): break
        else:
            read=True
            random_seed = numpy.random.random_integers(2024, size=8)
            res = function.wl_mdos(meanE,deltae,nbin,random_seed,read)

 
        
        if (res!=0):
            print "ERROR: some error in the calculation (python message)"
            return
        
        if (i==0):
            xdata = function.logdosf[:,0].copy()
            data  = numpy.zeros( len(xdata) )
            data2 = data.copy()
        
        if (len(function.logdosf[:,0]) != len(xdata) ):
            
            print "error:",len(xdata),len(function.logdosf[:,0])
            return

        res = utilities.HistoFromLog(function.logdosf[:,1],function.delta)
        data += res
        data2 += res*res
    #
    t2 = time()

    print "Elapsed time for calculation:",t2-t1,"seconds"

    data = data/numpy.float(nsim)
    data2 = data2/numpy.float(nsim)
    #
    sigma = numpy.sqrt(data2 - data*data)
    sigma = sigma/numpy.sqrt( numpy.float(nsim) )
    #
    tmp = numpy.column_stack( ( xdata , data, sigma ) )
    numpy.savetxt( datadir+'data0.out', tmp, fmt='%.18e' )
    
    Conf_numerical_plot()

    return


# calculation for differents L

def numerical( h0, h, deltae, nbin, nsim, threshold ):
    # ... stampare tempo impiegato ...
    #
    function = ising1D.wl
    ising1D.system.h = h
    ising1D.system.h0 = h0
    function.threshold = threshold
    #
    for L in range(20,121,20):
        ising1D.system.h = L
        meanE = system.E0(h0,h,L)
        #
        for i in range(nsim):
            j = 0
            if (i==0):
                read=False
                # repeat calculation if i'm not able to find an initial conf
                while (j<10):
                    j = j+1
                    random_seed = numpy.random.random_integers(2024, size=8)
                    res = function.wl_mdos(meanE,deltae,nbin,random_seed,read)
                    if (res==0): break
            else: 
                read=True
                random_seed = numpy.random.random_integers(2024, size=8)
                res = function.wl_mdos(meanE,deltae,nbin,random_seed,read)
        
            if (res!=0):
                print "ERROR: some error in the calculation (python message)"
                return
        
            if (i==0):
                xdata = function.logdosf[:,0].copy()
                data  = numpy.zeros( len(xdata) )
                data2 = data.copy()
                
            if (len(function.logdosf[:,0]) != len(xdata) ):
                    print "error:",len(xdata),len(function.logdosf[:,0])
                    return

            res = utilities.HistoFromLog(function.logdosf[:,1],function.delta)
            data += res
            data2 += res*res

    data = data/numpy.float(nsim)
    data2 = data2/numpy.float(nsim)
    #
    sigma = numpy.sqrt(data2 - data*data)
    sigma = sigma/numpy.sqrt( numpy.float(nsim) )
    #
    tmp = numpy.column_stack( ( xdata , data, sigma ) )
    numpy.savetxt( datadir+'data0.out', tmp, fmt='%.18e' )

    Conf_numerical_plot()

    return


def Conf_numerical_plot():
    data = numpy.loadtxt( datadir+'data0.out')

    fig = pp.figure()
    ax = fig.add_subplot(111)

    ax.plot(data[:,0],data[:,1],'-')

    ax.errorbar(data[:,0],data[:,1],yerr = 3.0*data[:,2],fmt=None)

    return
