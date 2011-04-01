
import sys
sys.path.append('../libs/')

import system
import numpy

import matplotlib.pyplot as pp
from matplotlib import rc

rc('text', usetex=True)


########################################
#### plot of single particle energy ####
########################################
def PlotEpsilon( h = [ 0.5 , 1.0 , 1.5 ] ):
    j = 1.0
    L = 100
    k = system.kpoints(L,0,"all")
    lines = []
    ll = []
#    pp.title(r"$$")
    for i in range(len(h)):
        ek = system.epsilon(k, h[i] ,j )
        lines.append( pp.plot(k,ek) )
        ll.append( "".join( ['$ h =',`h[i]`,'$' ] ) )
    pp.xlabel('$k$', fontsize=16)
    pp.ylabel(r'$\epsilon_k/J$',fontsize=16)
    pp.legend( lines, ll, loc=0)
    pp.xlim( -numpy.pi,numpy.pi )
    pp.show()
    return


########################################
#### plot of single particle energy ####
########################################
def PlotSigmaz( h = [ 0.5 , 1.0 ,1.5 ] ):
    import ising1D

    L = 100
    k = system.kpoints(L,0,"all")
    lines = []
    ll = []
    ising1D.system.l = L

    for i in range(len(h)):
        ising1D.system.h = h[i]
        ek = []
        for j in range( len(k) ):
            ek.append( ising1D.system.sigmaz( k[j] ) )

        lines.append( pp.plot(k,ek) )
        ll.append( "".join( ['$ h =',`h[i]`,'$' ] ) )

    pp.xlabel('$k$', fontsize=16)
    pp.ylabel(r'$m_k$',fontsize=16)
    pp.legend( lines, ll, loc=0)
    pp.xlim( -numpy.pi,numpy.pi )
    pp.show()
    return

##########################################
### trial to calculate m max and min   ###
### at given energy                    ###
##########################################
def PlotMaxMag( ):
    import ising1D
    
    h = 0.5
    L = 500

    k = system.kpoints(L,0,">")

    ising1D.system.L = L
    ising1D.system.h = h
    
    egs_p = ising1D.system.egs()
    egs_m = egs_p
    mgs_p = ising1D.system.mgs()
    mgs_m = mgs_p
    
    edata_p = [ egs_p ]
    edata_m = [ egs_m ]
    mdata_p = [ mgs_p ]
    mdata_m = [ mgs_m ]

    for i in range(len(k)):
        
        egs_p += 2.0*ising1D.system.energy( k[i], h, 1.0 )/numpy.float(L)
        mgs_p += 2.0*ising1D.system.sigmaz( k[i] )/numpy.float(L)

        egs_m += 2.0*ising1D.system.energy( k[len(k)-1-i], h, 1.0 )/numpy.float(L)
        mgs_m += 2.0*ising1D.system.sigmaz( k[len(k)-1-i] )/numpy.float(L)

        edata_p.append( egs_p )
        edata_m.append( egs_m )

        mdata_p.append( mgs_p )
        mdata_m.append( mgs_m )

    pp.plot(edata_p,mdata_p)
    pp.plot(edata_m,mdata_m)
    pp.show()
    
    return


#########################################
#### plot of energy after the quench ####
#########################################
def PlotE0():
    h0 = 0.05*numpy.arange(40)
    h  = [0.5,1.0,1.5]
    L  =  100
    lines = []
    ll = []
    pp.title( "".join( ["$L =",`L`,"$"] ),fontsize=20 )
    for i in range(len(h)):
        en = []
        for j in range(len(h0)):
            en.append( system.E0(h0[j],h[i],L) )
        lines.append(pp.plot(h0,en))
        ll.append( "".join( ["$h =",`h[i]`,"$"] ) )

    lines.append(pp.plot(h0,system.Egs(h0,L),'r--'))
    ll.append("$E_{\mathrm{gs}}/L$")
    pp.xlabel("$h^0/J$",fontsize=16)
    pp.ylabel("$E_0/(J L)$",fontsize=16)

    pp.annotate("",(0.5,system.E0(0.5,1.5,L)) ,(0.5,system.Egs(0.5,L)), arrowprops=dict(arrowstyle="->") )

    pp.annotate("",(1.5,system.E0(1.5,0.5,L)) ,(1.5,system.Egs(1.5,L)), arrowprops=dict(arrowstyle="->") )

    pp.legend( lines , ll, loc=0)
    pp.show()
    return

###########################################
#### plot of excess energy after the   ####
#### quench                            ####
###########################################
def PlotExEnergy():
    h0 = 0.05*numpy.arange(40)
    h = [0.5,1.0,1.5]
    L = 100

    lines = []
    ll = []
    pp.title( "".join( ["$L =",`L`,"$"] ),fontsize=20 )

    for i in range(len(h)):

        egstate = system.Egs( h[i], L )
        en = []
        for j in range(len(h0)):
            en.append( system.E0(h0[j],h[i],L) - egstate )
        lines.append(pp.plot(h0,en))
        ll.append( "".join( ["$h =",`h[i]`,"$"] ) )

    pp.xlabel("$h^0/J$",fontsize=16)
    pp.ylabel(r"$(E_0-E_{\mathrm{gs}})/(J L)$",fontsize=16)

    pp.legend( lines , ll, loc=0)
    pp.show()

    return



def PlotCoef():
    h0 = [0.5 , 1.0 , 1.5]
#attenzione che hai copiato solo il link
    h = h0
    L = 100
    lines = []
    ll = []
    k = system.kpoints(L,0,">")
    for i in range(len(h0)):
        for j in range(len(h)):
            if (h0[i] < h[j]):
                c = system.Coef( h0[i] , h[j] , L )
                lines.append( pp.plot( k,c ) )
                ll.append( "".join(["$h^0 =",`h0[i]`,"$, $h =",`h[j]`,"$" ]) )

    pp.title(''.join(["$L = ",`L`,"$"]))
    pp.xlabel("$k$",fontsize=16)
    pp.ylabel("$z_k$",fontsize=16)
    pp.legend( lines, ll, loc=0 )
    pp.yscale('log')
    pp.xlim(0,numpy.pi)
    pp.show()
    return


def PlotProj():
    h = [ 0.5, 1.0, 1.5 ]
    h0max = 2.0
    h0min = 0.0
    npoints = 50
    L = 100
    h0 = h0min + numpy.arange(npoints)*(h0max-h0min)/numpy.float(npoints)

    lines = []
    ll = []

    for i in range(len(h)):
        xdata = []
        ydata = []
        for j in range(npoints):
            xdata.append( h0[j] )
            ydata.append( system.CoefGS( h0[j], h[i], L ) )
        lines.append(pp.plot(xdata,ydata))
        ll.append( "".join(["$h = ",`h[i]`,"$"] ))

    pp.legend( lines, ll, loc=0 )

    pp.title(''.join(["$L = ",`L`,"$"]))
    pp.xlabel('$h^0$')
    pp.ylabel(r'$ | \left< \Psi | \Psi^0 \right> |^2$')
    pp.show()

    return



def trial(h,l):
    print system.Egs( h, l )
