
import numpy
import math

def selSubspace( data , parityArray, paritySel ):
    """Select a subspace of the data with a given parity"""
    out = []
    for i in range( len(parityArray) ):
        if ( parityArray[i] == paritySel):
            out.append( data[i,:] )

    return numpy.array( out )


def HistoFromLog( logdos, delta ):
    """Create normalized histo starting from log of dos"""
    beta = 0.0
    histo = numpy.zeros( len(logdos) )
    
    for i in range( len(logdos) ):
        histo[i] = math.exp( logdos[i] )
        beta += histo[i]

    histo = histo / (beta*delta)

    return numpy.array( histo )

def HistoFromData ( data, xmin, xmax, deltax ):
    """Create an histo using data"""
    n = int((xmax - xmin)/deltax)
    histo = numpy.zeros( n , dtype = numpy.int )
    xdata = xmin + deltax*( numpy.arange(n)+0.5 )
    
    for i in range(len(data)):
        ibin = int( (data[i] - xmin) / deltax )
        histo[ ibin ] += 1

    histo_out = []
    xdata = []
    
    for i in range(n):
        if (histo[ i ] != 0):
            xdata.append( xmin + deltax * ( numpy.float(i) + 0.5 ) )
            histo_out.append( numpy.float(histo[ i ])/ (numpy.float( len(data) )* deltax) )

    return numpy.array(xdata), numpy.array(histo_out)
    

def MeanValue ( data , weight ):
    """Calculate mean value and standard deviation using data and their weight"""
    mean = 0.0
    mean2 = 0.0
    ndata = len(data)
    
    for i in range(ndata):
        mean = mean + weight[i]*data[i]
        mean2 = mean2 + weight[i]*data[i]*data[i]
        
    sigma = numpy.sqrt((mean2-mean*mean))

    return mean,sigma

def ProbDist ( obs , calpha , nbin ):
    omin = obs.min()
    omax = obs.max()
    deltao = (omax-omin)/numpy.float(nbin)
    p = numpy.zeros( nbin )
    odata = omin+(numpy.arange(nbin)+0.5)*deltao
    for i in range(len(calpha)):
        ibin = numpy.int( (obs[ i ] - omin) / deltao )
        if (ibin==nbin):
            ibin = nbin -1
        p[ ibin ] = p[ ibin ] + calpha[ i ]
        
    p = p / deltao

    return numpy.array(odata),numpy.array(p)


def select( data, column, xmin, xmax ):
    """Select data of a specific column with value inside the interval [xmin,xmax]

    Usage:
          select( data, column, xmin, xmax)
    Where:
          data        --> array
          column      --> column to select
          xmin,xmax   --> min and max values of the interval
    """
    out = []
    for i in range( numpy.shape( data )[0] ):
        if (( data[i, column] >= xmin ) and ( data[i,column] <= xmax )):
            out.append(data[i,:])
    return numpy.array( out )


def sum( data, xColumn, dosColumn, fluctColumn = -1 ):
    """Calculate the marginal dos, optionally it calculate also the fluctuations

    Usage:
          sum( data, xColumn, dosColumn, fluctColumn)
    Where:
          data        --> 2D array whear each column correspond to one coordinate
          xColumn     --> coordinate to mantain in the marginal calculation
          dosColumn   --> 1D array with log of DOS
          flucColumn  --> if grater than 0 calculates also the fluctations
    remember that the result is \log ( \sum \exp( \log DOS ) )"""
    xval = []
    sum = []
    xval.append( data[ 0 , xColumn ] )
    sum.append( 0.0e0 )
    shift = data[ 0, dosColumn ]
    if (fluctColumn >= 0):
        fluctSum = []
        fluctSum.append( 0.0e0 )

    for i in range( numpy.shape(data)[0] ):
        j = 0
        while ( xval[j] != data[ i, xColumn ] ):
            j += 1
            if (j > ( len(xval) -1 ) ):
                xval.append( data[ i, xColumn ] )
                sum.append( 0.0e0 )
                if ( fluctColumn >= 0):
                    fluctSum.append( 0.0e0 )
                    
        sum[ j ] += math.exp( data[ i, dosColumn ] - shift )
        if ( fluctColumn >= 0):
            fluctSum[ j ] += (math.exp( data[ i, dosColumn ] - shift )* data[ i, fluctColumn ] )**2

    shift = math.log( sum[ 0 ] )

    out = []
    if (fluctColumn < 0):
        for i in range( len( sum ) ):
            out.append( [xval[i], math.log( sum[i] ) - shift] )
    else:
        for i in range( len( sum ) ):
            out.append( [xval[i], math.log( sum[i] ) - shift, math.sqrt( fluctSum[i]/(sum[i]**2) + fluctSum[0]/(sum[0]**2)  )] )

    return numpy.array( out )

