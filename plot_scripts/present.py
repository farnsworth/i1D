
import sys

sys.path.append('../libs/')
sys.path.append('../')

from all import datadir

import numpy
import ising1D
import utilities
import system
import matplotlib.pyplot as pp
import os.path
from time import sleep

from matplotlib import rc
rc('text', usetex=True)


def video( h , L , emin , emax , nbin ):
    #
    ising1D.system.l = L
    ising1D.system.h = h
    ising1D.wl.threshold = 0.2
    random_seed = numpy.random.random_integers(2024, size=8)
    #
    ising1D.wl.wl_edos( emin, emax, nbin, random_seed, False )
    #
    return

def produce_video():
    #
    #
    step = 0
    fig = pp.figure( figsize=[11.0,6.0])
    #
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.set_xlabel(r'$E/J$',fontsize=18 )
    ax1.set_ylabel(r'$\log \rho + \alpha $',fontsize=18)
    ax1.set_xlim(-1.66,1.66)
    ax1.set_xticks([-1.5,-1.0,-0.5,0,0.5,1.0,1.5])

    ax2.set_xlabel(r'$E/J$',fontsize=18)
    ax2.set_ylabel(r'$histogram$',fontsize=18)
#    ax2.set_ylim(0.0,1.5)
    ax2.axhspan(0.8, 1.2, facecolor='b', alpha=0.2)
    ax2.set_xlim(-1.66,1.66)
    ax2.set_xticks([-1.5,-1.0,-0.5,0,0.5,1.0,1.5])
    #
    frame = 0
    #
    f = numpy.exp(2.0)
    #
    while (True):

        step = step + 1
        i_swep = -20

        f = numpy.sqrt(f)
        if (step>1):
            del ax1.lines[:]

        if (step>10):break
        while (True):

            i_swep = i_swep + 20
            file_name = `step`+"-"+`i_swep`+".out"
            ex = os.path.isfile(file_name)
            

            if (ex):
                print file_name
                fig.suptitle( "step "+`step`+", loop "+`i_swep`+", log(f) = %.2e " %numpy.log(f),fontsize="20" )
#                fig.text(0.4, 0.95, "step "+`step`+", loop "+`i_swep`, size=20)
                if ((i_swep == 0)and(step == 1)):
                    file_name = "1-200.out"

                data = numpy.loadtxt(file_name)

                if (i_swep == 0)and(step==1):
                    data[:,1]=1
                    data[:,2]=0

                totdata = numpy.float( data[:,2].sum() )
                nbins = numpy.float(len(data[:,2]))
                if (i_swep != 0)or(step != 1):
                    mh = totdata/nbins
                else:
                    mh = 1

                histo = data[:,2]/mh
                ax1.plot(data[:,0],data[:,1],"b-")
                if (i_swep==0):
                    histo = histo*0
                ax2.plot(data[:,0],histo,"k-",linewidth=2.0 )
                ymax = numpy.int(data[:,1].max()/20.0)+1
                ax1.set_ylim(0,ymax*20.0 )
                ymax2 = numpy.int(histo.max())+1
                ymax2 = numpy.max( [ymax2,2]  )
                ax2.set_ylim(0,ymax2)
                pp.draw()
#
                frame += 1
                filename = str('%05d' % frame) + '.png'
                pp.savefig(filename, dpi=100)
                
                if ((i_swep==0)):
                    npause = 10
                    if (step<4):npause = 50
                    if (step==1):npause = 10
                    for i in range(npause):
                        frame += 1
                        filename = str('%05d' % frame) + '.png'
                        pp.savefig(filename, dpi=100)


#                sleep(0.5)
#                pp.draw()
            else:
                print "non trovato",file_name
                fig.suptitle( "step "+`step`+", loop "+`i_swep-20`+", log(f) = %.2e " %numpy.log(f),fontsize="20" )
                ax2.plot(data[:,0],histo,"k-",linewidth=2.0 )
                npause = 10
                if (step<3): npause = 50
                for i in range(npause):
                    frame += 1
                    filename = str('%05d' % frame) + '.png'
                    pp.savefig(filename, dpi=100)
                del ax2.lines[:]
                del fig.texts[-1]
                break

            del ax2.lines[:]
            del fig.texts[-1]

    os.system('ffmpeg -r 25 -i %05d.png out.avi')

    return
