
import sys
sys.path.append('../libs/')
sys.path.append('../')

import numpy
import math
import random as rnd
import ising1D
import utilities
import system
import matplotlib.pyplot as pp

from matplotlib import rc
rc('text', usetex=True)

###########################################
### attenzione che qualche riferimento  ###
### puo' essere andato a farsi benedire.###
###########################################
def comparison( emin, emax, deltam ):
    """Made a comparison between exact logdos and wl logdos.
    
    Usage:
          comparison(emin, emax, deltam)

    Where:
          emin,emax       --> minimum and maximum energy
          deltam          --> histogram interval"""


    ising1D.system.l = 12
    nbin = 11
    
#    seed = [ rnd.randrange(1,10000,1), rnd.randrange(1,10000,1) ]
#    wl.wl_mdos( emin, emax, deltam, seed )

    ising1D.exact.spectrum_calc()
    exactdata = utilities.selSubspace( ising1D.exact.spectrum, ising1D.exact.par, 0 )
    exactdata = utilities.select( exactdata, 0, emin, emax )

    seed = [ rnd.randrange(1,10000,1), rnd.randrange(1,10000,1) ]
    (mmin, mmax) = ising1D.wl.wl_mdos( emin, emax, deltam, seed )


    xdata, histo = utilities.HistoFromData( exactdata[:,1], mmin, mmax, deltam )

    histo = numpy.log(histo)
    histo = histo - histo[ 0 ]

    pp.plot( xdata, histo, 'ro')

    pp.plot( ising1D.wl.logdosf[:,0] , ising1D.wl.logdosf[:,1] ,'bo')
    pp.show()

    return


########################################################
### mean energy and its fluctuation in function of L ###
### for a quench                                     ###
########################################################
def fluct( ):
    """Mean energy and its fluctuation in function of L"""

    lmin = 6
    lmax = 51
    ldelta = 4
    data = []

    h = [ 0.5, 1.0, 1.5 ]
    h0 = [ 0.5, 1.0, 1.5 ]

    pp.yscale('log')
    pp.xscale('log')

    lines= []
    ll = []

    for i1 in range(len(h0)):
        for i2 in range(len(h)):
            if (h0[i1] != h[i2]):
                ising1D.system.h0 = h0[i1]
                ising1D.system.h = h[i2]
                data = []
    
                for l in range(lmin,lmax,ldelta):
                    ising1D.system.l = l
                    ising1D.quench.dist_calc(5)
                    a = ising1D.quench.mean_en
                    b = ising1D.quench.sigma_en
                    data.append( [l, numpy.float(a), numpy.float(b) ] )
    #                print ising1D.quench.mean_en, system.E0(h0,h, l)

                data = numpy.array(data)

                fit = numpy.polyfit(numpy.log(data[2:,0]),numpy.log(data[2:,2]),1)
                print fit
                xdata = numpy.arange(3,60,1)
                ydata = numpy.exp(fit[0]*numpy.log(xdata)+fit[1])

                lines.append( pp.plot( data[:,0], data[:,2],'o' ) )
#                pp.plot( data[:,0], data[:,2],'o' )
                pp.plot( xdata, ydata, 'b-')
                ll.append(r"".join(["$h^0 = ",`h0[i1]`,"\quad h = ",`h[i2]`,r" , \quad \alpha = ", "%.3f" %fit[0],"$"] ))
    
    pp.legend( lines, ll, loc = 0 )
    pp.xlabel('$L$')
    pp.ylabel(r'$\frac{\sqrt{<E^2>_D - <E>_D^2>}}{L}$')
#    pp.title( "".join( ["$h_0 = ",`h0`,"\qquad h =",`h`,"$"] ) )
    pp.show()

    return data


#################################################################
### plot for some L the c_alpha values in function of energy  ###
#################################################################
def PlotCalphaSeries(h0,h):
    """Mean energy and its fluctuation in function of L"""
    from matplotlib.ticker import MaxNLocator

# number of plots
    n = 2
# starting length
    lmin = 8
# steps of length
    ldelta = 8

    ising1D.system.h0 = h0
    ising1D.system.h = h

# f1 is the plot of c_\alpha
# f2 is the plot of the associated distribution
    f1,axarr1 = pp.subplots(n+1, sharex=True )
#    f2,axarr2 = pp.subplots(n+1, sharex=True, sharey=True )

#delete y ticks for the plots
    f1.subplots_adjust(hspace=0)
    pp.setp([a.get_xticklabels() for a in  axarr1[:n]],Visible=False)

#    f2.subplots_adjust(hspace=0)
#    pp.setp([a.get_xticklabels() for a in  axarr2[:n]],Visible=False)


    for i in range(0,n+1):
        l = lmin + i*ldelta
        print l
        ising1D.system.l = l
        ising1D.quench.calpha2_calc()
        c2 = ising1D.quench.calpha2
        meanE = system.E0(h0,h,l)

        axarr1[i].yaxis.set_major_locator(MaxNLocator(3))
        axarr1[i].plot( c2[:,0],c2[:,1],'o',label='L = '+str(l) )
        axarr1[i].legend()
        axarr1[i].axvline( meanE, linewidth = 2)
        
#        data = utilities.ProbDist( c2[:,0] , c2[:,1] , 50 )
#        axarr2[i].plot( data[0], data[1], 'o',label='L = '+str(l) )
#        axarr2[i].legend()

    pp.suptitle( ''.join( ['$h^0 = ',`h0`, ' \quad h = ',`h`,'$' ] ) )
    pp.xlabel('$E_a$')
    pp.ylabel(r'$|c_\alpha|^2$')
    pp.show()

    return

################################################################
# plot a single distribution of c_\alpha in a way to use uge L #
################################################################
# da cambiare #
###############
def PlotCalpha(h0,h,l):
    ising1D.system.h0 = h0
    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.quench.calpha2_calc()
    meanE = system.E0(h0,h,l)
    c2 = ising1D.quench.calpha2

    pp.plot( c2[:,0],c2[:,1],'o')    
    pp.xlabel(r"$E_\alpha$")
    pp.ylabel(r"$|c_\alpha|^2$")
    pp.axvline( meanE, linewidth=2 )

    pp.show()

    return

#####################################################################
# plot distribution generated by calpha without a uge use of memory #
#####################################################################
def PlotDist(h0,h,l,nbin):
    ising1D.system.h0 = h0
    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.quench.dist_calc(nbin)
    meanE = system.E0(h0,h,l)
    print meanE, ising1D.quench.mean_en
    dist = ising1D.quench.dist[:]

    lines = []
    ll = []

    lines.append(pp.plot( dist[:,0],dist[:,1],'bo-'))
    
    ll.append("".join(["$L = ",`l`,"$"] ) )

    pp.xlabel(r"$E_\alpha$")
    pp.ylabel(r"$P_\mathrm{D}$")
    pp.axvline( meanE, linewidth=2, color='r' )
    pp.axvline( ising1D.quench.mean_en, linewidth=2, color='b')
    pp.title(r"".join(['$h^0 = ',`h0`,'\qquad h = ',`h`,'$']))
    pp.legend( lines, ll, loc = 0 )
    pp.show()

    return


##############################################################
# plot the thermal mag canonical avarage in function of      #
# the mean energy canonical average                          # 
##############################################################
def PlotCanonicalMag(h , l):

    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA

    params = {'axes.labelsize': 25,
          'text.fontsize': 10,
          'legend.fontsize': 20,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15}
    pp.rcParams.update(params)


    ising1D.system.h = h
    ising1D.system.l = l

    h0max = 4.0
    h0min = 0.0

    gsen = ising1D.system.egs()

    # because the for start from 1 and we need 51 points
    npoints = 52

    # -h is the maximum obtainable energy with a quench
    # instead 0 is the maximum energy with thermal exitations
    emax = max( -h, -1.0 )
    delta = (emax - gsen )/numpy.float(npoints)
    en = []
    temp = []
    mag = []
    
    for i in range(1,npoints):
        en.append(gsen + delta*numpy.float(i) )
        temp.append( ising1D.canonical.find_temperature( en[i-1] , 0.000001 ) )
        mag.append( ising1D.canonical.thermal_magnetization( temp[i-1] ) )


    en2 = []
    mag2 = []
    delta = (h0max - h0min)/numpy.float(npoints)
    
    for i in range(0,npoints):
        h0 = delta*i + h0min
        ising1D.system.h0 = h0
        en2.append( system.E0(h0,h,l) )
        mag2.append( ising1D.quench.long_time_mag( ) )
        

    ax = host_subplot(111,axes_class=AA.Axes )

    lines = []
    lines.append(ax.plot(en, mag, linewidth=2.0 ) )
    lines.append(ax.plot(en2,mag2,"ro", linewidth=2.0) )
    ll = [r"$\langle M_z \rangle_\mathrm{c}$",r"$\langle M_z \rangle_\mathrm{D}$"]
#    ax.set_yticks([0.2,0.3,0.4,0.5])
    ax.set_yticks( [0.5,0.55,0.6] )


    ax.set_xlabel(r'$E_0$')
    ax.set_ylabel(r'$M_z$')

#    [i.set_fontsize(20) for i in ax.get_xticklabels() + ax.get_yticklabels()]


    # attenzione la temperatura che stampi non e' quella giusta
    # e' la meta' o il doppio di quella esatta, il fatto e' che dentro
    # la funzione di Fermi mandi temperatura
    # a way to have a second axes with different ticks
    ax2 = ax.twin() # ax2 is responsible for "top" axis and "right" axis
    tmax = round(temp[npoints-2],1)
    tmin = 0.0
    delta = round((tmax-tmin)/6.0,1)
    tmax = tmax + delta
    ax2.set_xticks( [  (ising1D.system.egs() + ising1D.canonical.thermal_ex_energy(round(tmin + delta*numpy.float(i),2) ))   for i in range(1,8) ] )
    ax2.set_xticklabels([r"$%.1f$" %(round(tmin + delta*numpy.float(i),2) ) for i in range(1,8) ] )

    ax2.set_xlabel("$T_{eff}$",fontsize=30)
    ax2.axis["right"].major_ticklabels.set_visible(False)

#    ax.annotate( ('L = '+`l`+", h ="+`h`) , xy=(.8, .8),  xycoords='axes fraction',
#                horizontalalignment='center', verticalalignment='center')
 
    ax.set_xlim( xmin = gsen, xmax = emax )
#    ax.set_xlim( xmin = gsen, xmax = -0.65 )
    leg = ax.legend( lines, ll, loc = 0 )

    for t in leg.get_texts():
        t.set_fontsize(20)    # the legend text fontsize


    pp.draw()
    pp.show()





##############################################################
# plot the thermal kinks canonical avarage in function of    #
# the mean energy canonical average                          # 
##############################################################
def PlotCanonicalKinks(h , l):

    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA


    ising1D.system.h = h
    ising1D.system.l = l

    h0max = 4.0
    h0min = 0.0

    gsen = ising1D.system.egs()

    # because the for start from 1 and we need 51 points
    npoints = 52

    # -h is the maximum obtainable energy with a quench
    # instead 0 is the maximum energy with thermal exitations
    emax = max( -h, -1.0 )
    delta = (emax - gsen )/numpy.float(npoints)
    en = []
    temp = []
    kinks = []
    
    for i in range(1,npoints):
        en.append(gsen + delta*numpy.float(i) )
        temp.append( ising1D.canonical.find_temperature( en[i-1] , 0.000001 ) )
        kinks.append( ising1D.canonical.thermal_kinks( temp[i-1] ) )


    en2 = []
    kinks2 = []
    delta = (h0max - h0min)/numpy.float(npoints)
    
    for i in range(0,npoints):
        h0 = delta*i + h0min
        ising1D.system.h0 = h0
        en2.append( system.E0(h0,h,l) )
        kinks2.append( ising1D.quench.long_time_kinks( ) )

    ax = host_subplot(111,axes_class=AA.Axes )

    lines = []
    lines.append(ax.plot(en, kinks ))
    lines.append(ax.plot(en2, kinks2, "ro-"))
    ll = [r"$\left< N \right>_c (T_\mathrm{eff})$",r"$\left< N \right>_t$"]

    # a way to have a second axes with different ticks
    ax2 = ax.twin() # ax2 is responsible for "top" axis and "right" axis
    tmax = round(temp[npoints-2],1)
    tmin = 0.0
    delta = round((tmax-tmin)/6.0,1)
    tmax = tmax + delta
    ax2.set_xticks( [  (ising1D.system.egs() + ising1D.canonical.thermal_ex_energy(round(tmin + delta*numpy.float(i),2) ))   for i in range(0,8) ] )
    ax2.set_xticklabels(["$%.1f$" %(round(tmin + delta*numpy.float(i),2) ) for i in range(0,8) ] )

    ax2.set_xlabel("$T_{eff}$")
    ax2.axis["right"].major_ticklabels.set_visible(False)

    ax.annotate( ('L = '+`l`+", h ="+`h`) , xy=(.2, .8),  xycoords='axes fraction',
                horizontalalignment='center', verticalalignment='center')

    ax.set_xlim( xmin = gsen, xmax = emax )
    pp.legend( lines, ll, loc = 0 )
    pp.xlabel('$E_0$')
    pp.ylabel('$N$')

    pp.draw()
    pp.show()


#############################################################
# plot a single distribution of c_\alpha in function of mag #
#############################################################
def PlotMagalpha(h0,h,l):
    ising1D.system.h0 = h0
    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.quench.calpha2_calc()

    ising1D.quench.magalpha_calc()
    obs = ising1D.quench.obs.copy()

    ising1D.quench.ealpha_calc()
    ealpha = ising1D.quench.obs.copy()


    equench = system.E0(h0,h, l)
    temp = ising1D.canonical.find_temperature( equench , 0.000001 )
    canmag = ising1D.canonical.thermal_magnetization( temp )

    timemag = ising1D.quench.long_time_mag()
    
    c2 = ising1D.quench.calpha2


    pp.plot( obs ,c2,'o')    
    pp.xlabel(r"${M_z{}}_{\alpha \alpha}$")
    pp.ylabel(r"$|c_\alpha|^2$")
    pp.axvline( canmag, linewidth=2 )
    pp.axvline( timemag, linewidth=2 )

    pp.show()

    return


#####################################################################
# plot distribution generated by calpha without a uge use of memory #
# in function of magnetization                                      # 
#####################################################################
def PlotMagDist(h0,h,l,nbin, deltaen ):
    ising1D.system.h0 = h0
    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.quench.mag_dist_calc(nbin)

    meanE = system.E0(h0,h,l)
    
    dist = ising1D.quench.dist.copy()

    ising1D.quench.sigma_mag_quench2( meanE , deltaen , nbin )
    dist2 = ising1D.quench.dist.copy()

    temp = ising1D.canonical.find_temperature( meanE , 0.000001 )
    canmag = ising1D.canonical.thermal_magnetization( temp )

    timemag = ising1D.quench.long_time_mag()

    print "Energy", meanE
    print "Can av",canmag
    print "Time av",timemag

    lines = []
    ll = []

    fig = pp.figure( )
    ax = fig.add_subplot(111)

    lines.append(ax.plot( dist[:,0],dist[:,1],'ro-'))
    lines.append( pp.axvline( timemag, linewidth=1.5, color = "r", linestyle="-") )
    lines.append(ax.plot( dist2[:,0],dist2[:,1],'bo-'))
    lines.append( ax.axvline( canmag, linewidth=1.5, linestyle="--" ) )

    
    ll.append(r"$P_{\mathrm{D}}$" )
    ll.append( r"$\left< M_z \right>_D$")
    ll.append(r"".join(["$P_{\mathrm{mc}}, \delta E/J =","%.2f" %deltaen," $"]))
    ll.append(r"$\left< M_z \right>_c (T_{eff}) $")

    ax.set_xlabel(r"$M_z$",fontsize=25)
    ax.set_ylabel(r"$P (M_z)$",fontsize=25)


 


#    pp.title(r"".join(['$h^0 = ',`h0`,'\qquad h = ',`h`,'\qquad L = ',`l`,'$']))
    leg = ax.legend( lines, ll, loc = 0 )

    [i.set_fontsize(20) for i in ax.get_xticklabels() + ax.get_yticklabels()]

    for t in leg.get_texts():
        t.set_fontsize(20)    # the legend text fontsize
        

    pp.show()

    return


def PlotKinksDist(h0, h, l, nbin, deltaen):

    ising1D.system.h0 = h0
    ising1D.system.h = h
    ising1D.system.l = l

    ising1D.quench.kinks_dist_calc(nbin)

    meanE = system.E0(h0,h,l)
    
    dist = ising1D.quench.dist.copy()

    ising1D.quench.sigma_kink_quench2( meanE , deltaen , nbin )
    dist2 = ising1D.quench.dist.copy()

    temp = ising1D.canonical.find_temperature( meanE , 0.000001 )
    cankinks = ising1D.canonical.thermal_kinks( temp )

    timekinks = ising1D.quench.long_time_kinks()

    print "Energy", meanE
    print "Can av",cankinks
    print "Time av",timekinks

    lines = []
    ll = []

    lines.append(pp.plot( dist[:,0],dist[:,1],'ro-'))
    lines.append(pp.plot( dist2[:,0],dist2[:,1],'bo-'))

    ll.append(r"$P_{\mathrm{D}}$" )
    ll.append(r"".join(["$P_{\mathrm{mc}}, \delta E/J =","%.3f" %deltaen," $"]))


    pp.xlabel(r"$N$")
    pp.ylabel(r"$P (N)$")
    lines.append( pp.axvline( cankinks, linewidth=1.5, linestyle="--" ) )
    lines.append( pp.axvline( timekinks, linewidth=1.5, color = "r", linestyle="-") )
    ll.append(r"$\left< N \right>_c (T_{eff}) $")
    ll.append( r"$\left< N \right>_D$")
    pp.title(r"".join(['$h^0 = ',`h0`,'\qquad h = ',`h`,'\qquad L = ',`l`,'$']))
    pp.legend( lines, ll, loc = 0 )
    pp.show()

    return



##########################################
### plot the dist for two differents l ###
##########################################
def PlotMagDist2(h0,h,l,nbin):
    ising1D.system.h0 = h0
    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.quench.mag_dist_calc(nbin)

    meanE = system.E0(h0,h,l)
    
    dist = ising1D.quench.dist.copy()

    temp = ising1D.canonical.find_temperature( meanE , 0.000001 )
    canmag = ising1D.canonical.thermal_magnetization( temp )

    timemag = ising1D.quench.long_time_mag()

    ising1D.system.l = 20
    ising1D.quench.mag_dist_calc( nbin )
    dist2 = ising1D.quench.dist.copy()


    print "Energy", meanE
    print "Can av",canmag
    print "Time av",timemag

    lines = []
    ll = []

    lines.append(pp.plot( dist[:,0],dist[:,1],'ro--'))
    ll.append("".join(["$L = ",`l`,"$"] ) )


    lines.append(pp.plot( dist2[:,0],dist2[:,1],'ro-'))
    ll.append("".join(["$L = 30$"] ) )

    pp.xlabel(r"$M_z$")
    pp.ylabel(r"$P_\mathrm{D}(M_z)$")
    lines.append( pp.axvline( canmag, linewidth=1.5, linestyle="--" ) )
    lines.append( pp.axvline( timemag, linewidth=1.5, color = "r", linestyle="-") )
    ll.append(r"$\left< M_z \right>_c (T_{eff}) $")
    ll.append( r"$\left< M_z \right>_D$")
    pp.title(r"".join(['$h^0 = ',`h0`,'\qquad h = ',`h`,'$']))
    pp.legend( lines, ll, loc = 0 )
    pp.show()

    return

def PlotKinksDist2(h0,h,l,nbin):
    ising1D.system.h0 = h0
    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.quench.kinks_dist_calc(nbin)

    meanE = system.E0(h0,h,l)
    
    dist = ising1D.quench.dist.copy()

    temp = ising1D.canonical.find_temperature( meanE , 0.000001 )
    cankinks = ising1D.canonical.thermal_kinks( temp )

    timekinks = ising1D.quench.long_time_kinks()

    ising1D.system.l = 20
    ising1D.quench.kinks_dist_calc( nbin )
    dist2 = ising1D.quench.dist.copy()


    print "Energy", meanE
    print "Can av",cankinks
    print "Time av",timekinks

    lines = []
    ll = []

    lines.append(pp.plot( dist[:,0],dist[:,1],'ro--'))
    ll.append("".join(["$L = ",`l`,"$"] ) )


    lines.append(pp.plot( dist2[:,0],dist2[:,1],'ro-'))
    ll.append("".join(["$L = 30$"] ) )

    pp.xlabel(r"$N$")
    pp.ylabel(r"$P_\mathrm{D}(N)$")
    lines.append( pp.axvline( cankinks, linewidth=1.5, linestyle="--" ) )
    lines.append( pp.axvline( timekinks, linewidth=1.5, color = "r", linestyle="-") )
    ll.append(r"$\left< N \right>_c (T_{eff}) $")
    ll.append( r"$\left< N \right>_D$")
    pp.title(r"".join(['$h^0 = ',`h0`,'\qquad h = ',`h`,'$']))
    pp.legend( lines, ll, loc = 0 )
    pp.show()

    return




def Plotmgs(l):

    ising1D.system.l = l
    npoints = 50
    hmax = 3.0
    hmin = 0.0
    deltah = (hmax-hmin)/numpy.float(npoints)
    harray = []
    marray = []

    for i in range(50):
        h = hmin + numpy.float(i) * deltah
        harray.append( h )
        ising1D.system.h = h
        marray.append( ising1D.system.mgs( ) )

    fig = pp.figure( figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.set_ylim([0,1])
    ax.plot(harray,marray, linewidth=2.5)
    ax.set_xlabel('$h/J$',fontsize=30)
    ax.set_ylabel('$M_z$',fontsize=30 )
    ax.set_position([0.12,0.18,0.8,0.7])

    ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0])
    ax.set_xticks([0,1,2,3])
    [i.set_fontsize(25) for i in ax.get_xticklabels() + ax.get_yticklabels()]
    [i.set_linewidth(1.5) for i in ax.spines.itervalues()]
    [i.set_markeredgewidth(1.5) for i in ax.get_xticklines() + ax.get_yticklines() ]

#    ax.set_title('Ground state transverse magnetization')

    return

def Plotkinksgs(l):

    ising1D.system.l = l
    npoints = 50
    hmax = 3.0
    hmin = 0.0
    deltah = (hmax-hmin)/numpy.float(npoints)
    harray = []
    marray = []

    for i in range(50):
        h = hmin + numpy.float(i) * deltah
        harray.append( h )
        ising1D.system.h = h
        marray.append( ising1D.system.kinksgs( ) )

    fig = pp.figure(figsize=(8,5))
    ax = fig.add_subplot(111)
    ax.plot(harray,marray,linewidth=2.5)
#    ax.set_aspect('auto',adjustable='box')
    ax.set_position([0.12,0.18,0.8,0.7])
    ax.set_xlabel('$h/J$',fontsize=30)
    ax.set_ylabel(r'$N_{\mathrm{d}}$',fontsize=30)
    ax.set_yticks([0,0.1,0.2,0.3,0.4,0.5])
    ax.set_xticks([0,1,2,3])
    [i.set_fontsize(25) for i in ax.get_xticklabels() + ax.get_yticklabels()]
    [i.set_linewidth(1.5) for i in ax.spines.itervalues()]
    [i.set_markeredgewidth(1.5) for i in ax.get_xticklines() + ax.get_yticklines() ]


#    ax.set_title('Ground state density of kinks')

    return


def PlotMagArray(h, l, ebin, mbin ):

    import matplotlib.cm as cm

    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.exact.mag_array_calc( ebin, mbin)
    
    data = ising1D.exact.array
    emax = ising1D.exact.emax_plot
    mmax = ising1D.exact.obsmax_plot

    fig = pp.figure()
    ax = fig.add_subplot(111)
    data = numpy.log(data+1)
# forse e' da cambiare questa riga
    cax = ax.imshow( numpy.transpose(data), origin="lower", interpolation="nearest", extent=(-emax,emax,-mmax,mmax),aspect='auto',cmap=cm.gray_r )
    cbar = fig.colorbar( cax )
    cbar.set_label(r"$\log ( n(E,M) + 1 ) $")
    ax.set_ylabel(r"$M_z$")
    ax.set_xlabel(r"$E/L$")
    ax.set_title( r"".join( ["$L = ",`l`,"$, $h = ",`h`,"$"] ) )



#    k = system.kpoints(l,0,">")
#    
#    egs_p = ising1D.system.egs()
#    egs_m = egs_p
#    mgs_p = ising1D.system.mgs()
#    mgs_m = mgs_p
#    
#    edata_p = [ egs_p ]
#    edata_m = [ egs_m ]
#    mdata_p = [ mgs_p ]
#    mdata_m = [ mgs_m ]
#
#    smaxlast = 1.5
#    sminlast = -1.5
#
#    for i in range(len(k)):
#        smaxk = -1.5
#        smink = 1.5
#        for j in range(len(k)):
#            sk = ising1D.system.sigmaz( k[j] )/numpy.float(l)
#            ek = ising1D.system.energy( k[j],h,1.0)/numpy.float(l)
#            if ( (sk < smaxlast)and(sk > smaxk) ):
#                smaxk = sk
#                emaxk = ek
#            if ( (sk > sminlast)and(sk < smink) ):
#                smink = sk
#                emink = ek
#
#        print smink, smaxk
#        egs_p += 2.0*emaxk
#        mgs_p += 2.0*smaxk
#        egs_m += 2.0*emink
#        mgs_m += 2.0*smink
#
#        edata_p.append( egs_p )
#        edata_m.append( egs_m )
#        mdata_p.append( mgs_p )
#        mdata_m.append( mgs_m )
#
#        smaxlast = smaxk
#        sminlast = smink
#
#    ax.plot(edata_p,mdata_p,'g--')
#    ax.plot(edata_m,mdata_m,'g--')
    ax.set_xlim(-emax,emax)
    ax.set_ylim(-mmax,mmax)
    
    return


def PlotKinksArray(h, l, ebin, kinksbin ):

    import matplotlib.cm as cm

    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.exact.kinks_array_calc( ebin, kinksbin)
    
    data = ising1D.exact.array
    emax = ising1D.exact.emax_plot
    kinksmax = ising1D.exact.obsmax_plot

    fig = pp.figure()
    ax = fig.add_subplot(111)
    data = numpy.log(data+1)
# forse e' da cambiare questa riga
    cax = ax.imshow( numpy.transpose(data), origin="lower", interpolation="nearest", extent=(-emax,emax,-kinksmax+0.5,kinksmax+0.5),aspect='auto',cmap=cm.gray_r )
    cbar = fig.colorbar( cax )
    cbar.set_label(r"$\log ( n(E,M) + 1 ) $")
    ax.set_ylabel(r"$N$")
    ax.set_xlabel(r"$E/L$")
    ax.set_title( r"".join( ["$L = ",`l`,"$, $h = ",`h`,"$"] ) )



def PlotSigmaMag(h0, h ):
    
    ising1D.system.h = h
    
    dearr = [0.1, 0.05, 0.01 ]

    fig = pp.figure()
    fig2 = pp.figure()

    ax = fig.add_subplot(111)
    ax2 = fig2.add_subplot(111)

    ax2.set_yscale('log')
    ax2.set_xscale('log')

    lines = []
    lines2 = []
    ll = []
    ll2 = []

#    alldata = []

    for deltae in dearr:
        xdata = []
        ydata = []
#mettere secondo numero a 61
        for i in range(20,61,2):
            xdata.append( i )
            ising1D.system.l = i
            energy = system.E0(h0,h,i)
            ydata.append( ising1D.quench.sigma_mag_quench2( energy, deltae , 1 ) )
        lines.append( ax.plot(xdata,ydata,'o-') )
        lines2.append( ax2.plot( xdata, ydata,'o'))
        ll.append( r"".join( ['$\delta E = $',"%.3f" %deltae] ) )

#### do a fit
        fit = numpy.polyfit(numpy.log(xdata[5:]),numpy.log(ydata[5:]),1)
        print fit
        xdata2 = numpy.arange(20,60,1)
        ydata2 = numpy.exp(fit[0]*numpy.log(xdata2)+fit[1])
        
        ax2.plot( xdata2, ydata2, '-')

        ll2.append( "".join( [r'$\delta E = ',' %.3f' %deltae, r', \alpha = %.3f' %fit[0], '$' ] ) )

        
    ax.set_title(r"".join( ["$h^0 = ",`h0`,", h = ",`h`, "$" ] ) )
    ax.set_xlabel(r"$L$")
    ax.set_ylabel(r"$\Delta M_z$")
    ax.legend( lines, ll, loc = 0 )

    ax2.set_title(r"".join( ["$h^0 = ",`h0`,", h = ",`h`, "$" ] ) )
    ax2.set_xlabel(r"$L$")
    ax2.set_ylabel(r"$\Delta M$")
    ax2.set_xlim( ax.get_xlim() )
    ax2.set_ylim( ax.get_ylim() ) 
#    ax2.set_xticks( ax.get_xticks() )
#    ax2.set_xticklabels( ax.get_xticklabels() )
    ax2.legend( lines2, ll2, loc = 0 )

    return

def PlotSigmaKinks(h0, h ):
    
    ising1D.system.h = h
    
    dearr = [0.1, 0.05, 0.01 ]

    fig = pp.figure()
    fig2 = pp.figure()

    ax = fig.add_subplot(111)
    ax2 = fig2.add_subplot(111)

    ax2.set_yscale('log')
    ax2.set_xscale('log')

    lines = []
    lines2 = []
    ll = []
    ll2 = []


    for deltae in dearr:
        xdata = []
        ydata = []
        for i in range(20,61,2):
            xdata.append( i )
            ising1D.system.l = i
            energy = system.E0(h0,h,i)
            ydata.append( ising1D.quench.sigma_kink_quench2( energy, deltae , 1 ) )
        lines.append( ax.plot(xdata,ydata,'o-') )
        lines2.append( ax2.plot( xdata, ydata,'o'))
        ll.append( r"".join( ['$\delta E = $',"%.3f" %deltae] ) )

#### do a fit
        fit = numpy.polyfit(numpy.log(xdata[5:]),numpy.log(ydata[5:]),1)
        print fit
        xdata2 = numpy.arange(20,60,1)
        ydata2 = numpy.exp(fit[0]*numpy.log(xdata2)+fit[1])
        
        ax2.plot( xdata2, ydata2, '-')

        ll2.append( "".join( [r'$\delta E = ',' %.3f' %deltae, r', \alpha = %.3f' %fit[0], '$' ] ) )

        
    ax.set_title(r"".join( ["$h^0 = ",`h0`,", h = ",`h`, "$" ] ) )
    ax.set_xlabel(r"$L$")
    ax.set_ylabel(r"$\Delta N$")
    ax.legend( lines, ll, loc = 0 )

    ax2.set_title(r"".join( ["$h^0 = ",`h0`,", h = ",`h`, "$" ] ) )
    ax2.set_xlabel(r"$L$")
    ax2.set_ylabel(r"$\Delta M$")
    ax2.set_xlim( ax.get_xlim() )
    ax2.set_ylim( 0.02,ax.get_ylim()[1] ) 
#    ax2.set_xticks( ax.get_xticks() )
#    ax2.set_xticklabels( ax.get_xticklabels() )
    ax2.legend( lines2, ll2, loc = 0 )

    return





def PlotMagArrayAllSpace(h, l, ebin, mbin ):

    import matplotlib.cm as cm

    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.exact.mag_array_calc_allspace( ebin, mbin)
#    ising1D.exact.mag_array_calc( ebin, mbin)
    
    data = ising1D.exact.real_array
    data2 = ising1D.exact.array
    emax = ising1D.exact.emax_plot
    mmax = ising1D.exact.obsmax_plot

#    tmp = data.min()
#    for i in range(ebin):
#        for j in range(mbin):
#            if ( data2[i,j] == 0 ):
#                data[i,j] = tmp
            
    data = numpy.log(data)

    fig = pp.figure()
    ax = fig.add_subplot(111)

    cax = ax.imshow( numpy.transpose(data), origin="lower", interpolation="nearest", extent=(-emax,emax,-mmax,mmax),aspect='auto',cmap=cm.gray_r )
#    cbar = fig.colorbar( cax )
#    cbar.set_label(r"$\log ( \rho(E,M) ) $")
    ax.set_ylabel(r"${M_z{}}_{\alpha \alpha}$",fontsize=20)
    ax.set_xlabel(r"$E_\alpha/(L J)$",fontsize=20)
    [i.set_fontsize(15) for i in ax.get_xticklabels() + ax.get_yticklabels()]
#    ax.set_title( r"".join( ["$L = ",`l`,"$, $h = ",`h`,"$"] ) )



    k = system.kpoints(l,0,">")
    
    egs_p = ising1D.system.egs()
    egs_m = egs_p
    mgs_p = ising1D.system.mgs()
    mgs_m = mgs_p
    
    edata_p = [ egs_p ]
    edata_m = [ egs_m ]
    mdata_p = [ mgs_p ]
    mdata_m = [ mgs_m ]

    for i in range(len(k)):
        
        egs_p += 2.0*ising1D.system.energy( k[i], h, 1.0 )/numpy.float(l)
        mgs_p += 2.0*ising1D.system.sigmaz( k[i] )/numpy.float(l)

        egs_m += 2.0*ising1D.system.energy( k[len(k)-1-i], h, 1.0 )/numpy.float(l)
        mgs_m += 2.0*ising1D.system.sigmaz( k[len(k)-1-i] )/numpy.float(l)

        edata_p.append( egs_p )
        edata_m.append( egs_m )

        mdata_p.append( mgs_p )
        mdata_m.append( mgs_m )

    ax.plot(edata_p,mdata_p,'g--')
    ax.plot(edata_m,mdata_m,'g--')
    ax.set_xlim(-emax,emax)
    ax.set_ylim(-mmax,mmax)
    
    return


def PlotMagArrayDifference(h, l, ebin, mbin ):

    import matplotlib.cm as cm

    ising1D.system.h = h
    ising1D.system.l = l
    ising1D.exact.mag_array_calc_allspace( ebin, mbin)

    data_all = ising1D.exact.array.copy()

    ising1D.system.l = 2*l
    ising1D.exact.mag_array_calc( ebin, mbin)
    data = ising1D.exact.array.copy()

    data = abs(numpy.log(data_all+1) - numpy.log(data+1))

    emax = ising1D.exact.emax_plot
    mmax = ising1D.exact.obsmax_plot

    fig = pp.figure()
    ax = fig.add_subplot(111)
    cax = ax.imshow( numpy.transpose(data), origin="lower", interpolation="nearest", extent=(-emax,emax,-mmax,mmax),aspect='auto',cmap=cm.gray_r )
    cbar = fig.colorbar( cax )
    cbar.set_label(r"$\delta $")
    ax.set_ylabel(r"$M_z$")
    ax.set_xlabel(r"$E/L$")
    ax.set_title( r"".join( ["$L = ",`l`,"$, $h = ",`h`,"$"] ) )

    return

#    k = system.kpoints(l,0,">")
#    
#    egs_p = ising1D.system.egs()
#    egs_m = egs_p
#    mgs_p = ising1D.system.mgs()
#    mgs_m = mgs_p
#    
#    edata_p = [ egs_p ]
#    edata_m = [ egs_m ]
#    mdata_p = [ mgs_p ]
#    mdata_m = [ mgs_m ]
#
#    for i in range(len(k)):
#        
#        egs_p += 2.0*ising1D.system.energy( k[i], h, 1.0 )/numpy.float(l)
#        mgs_p += 2.0*ising1D.system.sigmaz( k[i] )/numpy.float(l)
#
#        egs_m += 2.0*ising1D.system.energy( k[len(k)-1-i], h, 1.0 )/numpy.float(l)
#        mgs_m += 2.0*ising1D.system.sigmaz( k[len(k)-1-i] )/numpy.float(l)
#
#        edata_p.append( egs_p )
#        edata_m.append( egs_m )
#
#        mdata_p.append( mgs_p )
#        mdata_m.append( mgs_m )
#
#    ax.plot(edata_p,mdata_p,'g--')
#    ax.plot(edata_m,mdata_m,'g--')
#    ax.set_xlim(-emax,emax)
#    ax.set_ylim(-mmax,mmax)
#    
#    return







def PlotSigmaMagDiagonal( h0, h ):

    ising1D.system.h0 = h0
    ising1D.system.h = h


    data = []
    xdata = []
    
    for i in range(10,60,2):
        print "Simulation N = ",i
        ising1D.system.l = i
        xdata.append(i)
        ising1D.quench.mag_dist_calc(1)
        print "sigma",ising1D.quench.sigma_obs
        data.append( ising1D.quench.sigma_obs.copy() )

    fig = pp.figure()
    ax = fig.add_subplot(111)
    ax.plot(xdata, data)

    fig1 = pp.figure()
    ax1 = fig1.add_subplot(111)
    ax1.plot(xdata, numpy.log(data))

    fig2 = pp.figure()
    ax2 = fig2.add_subplot(111)
    ax2.plot(numpy.log(xdata), numpy.log(data))

    return


####
## verifica che i quench che partono da uno non termalizzano (forse)
###

def Trial( l):

    from mpl_toolkits.axes_grid1 import host_subplot
    import mpl_toolkits.axisartist as AA

    ising1D.system.h0 = 1.0
    ising1D.system.l = l

    hmax = 4.0
    hmin = 0.0

    gsen = ising1D.system.egs()

    # because the for start from 1 and we need 51 points
    npoints = 52

    en = []
    temp = []
    mag = []
    
#    for i in range(1,npoints):
#        en.append(gsen + delta*numpy.float(i) )
#        temp.append( ising1D.canonical.find_temperature( en[i-1] , 0.000001 ) )
#        mag.append( ising1D.canonical.thermal_magnetization( temp[i-1] ) )


    en2 = []
    mag2 = []
    delta = (hmax - hmin)/numpy.float(npoints)
    
    for i in range(0,npoints):
        h = delta*i + hmin
        ising1D.system.h = h
        en2.append( system.E0(1.0,h,l) )
        mag2.append( ising1D.quench.long_time_mag( ) )

######### qui devi mettere l'energia del quench
        en.append(  system.E0(1.0,h,l) )
        temp.append( ising1D.canonical.find_temperature( en[i] , 0.000001 ) )
        mag.append( ising1D.canonical.thermal_magnetization( temp[i-1] ) )
        

    ax = host_subplot(111,axes_class=AA.Axes )

    lines = []
    lines.append(ax.plot(en, mag, linewidth=2.0 ) )
    lines.append(ax.plot(en2,mag2,"ro", linewidth=2.0) )
    ll = [r"$\langle M_z \rangle_\mathrm{c}$",r"$\langle M_z \rangle_\mathrm{D}$"]


    ax.set_xlabel(r'$E$')
    ax.set_ylabel(r'$M_z$')

#    [i.set_fontsize(20) for i in ax.get_xticklabels() + ax.get_yticklabels()]


    # attenzione la temperatura che stampi non e' quella giusta
    # e' la meta' o il doppio di quella esatta, il fatto e' che dentro
    # la funzione di Fermi mandi temperatura
    # a way to have a second axes with different ticks
    ax2 = ax.twin() # ax2 is responsible for "top" axis and "right" axis
    tmax = round(temp[npoints-2],1)
    tmin = 0.0
    delta = round((tmax-tmin)/6.0,1)
    tmax = tmax + delta
    ax2.set_xticks( [  (ising1D.system.egs() + ising1D.canonical.thermal_ex_energy(round(tmin + delta*numpy.float(i),2) ))   for i in range(1,8) ] )
    ax2.set_xticklabels([r"$%.1f$" %(round(tmin + delta*numpy.float(i),2) ) for i in range(1,8) ] )

    ax2.set_xlabel("$T_{eff}$",fontsize=30)
    ax2.axis["right"].major_ticklabels.set_visible(False)

#    ax.annotate( ('L = '+`l`+", h ="+`h`) , xy=(.8, .8),  xycoords='axes fraction',
#                horizontalalignment='center', verticalalignment='center')
 
#    ax.set_xlim( xmin = gsen, xmax = emax )
#    ax.set_xlim( xmin = gsen, xmax = -0.65 )
    leg = ax.legend( lines, ll, loc = 0 )

    for t in leg.get_texts():
        t.set_fontsize(20)    # the legend text fontsize


    pp.draw()
    pp.show()
