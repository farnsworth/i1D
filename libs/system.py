
import numpy


def kpoints(L,parity,string):
    k = []
    delta = 2.0*numpy.pi/numpy.float(L)
#odd number of fermions
    if (parity != 0):
        if ( string == ">"):
            for i in range(1,L/2+1):
                k.append( numpy.float(i)*delta )
        else:
            for i in range(-L/2+1,L/2+1):
                k.append( numpy.float(i)*delta )
    else:
#even number of fermions
        delta0 = - numpy.pi + numpy.pi/numpy.float(L)
        if ( string == ">"):
            for i in range(L/2,L):
                k.append( delta0 + numpy.float(i)*delta )
        else:
            for i in range(0,L):
                k.append( delta0 + numpy.float(i)*delta )
    return numpy.array(k)


#single particle energy

def epsilon(k,h,j):
    temp = 2*numpy.sqrt(j*j + h*h -2.0*j*h*numpy.cos(k) )
    return temp

#energy after the quench
def E0(h0,h,L):
    j = 1.0
    k = kpoints(L,0,">")
    en = 0.0
    const = j*j+ h*h0

    for i in range(len(k)):
        if (epsilon(k[i],h0,j) == 0.0 ):
            print "warning, division by 0"
        else:
            en = en + (const - j*numpy.cos(k[i])*(h+h0))/epsilon(k[i],h0,j)
            
    en = -4.0*en/numpy.float(L)
    return en
    
#ground state energy
def Egs(h,L):
    j = 1.0
    k = kpoints(L,0,">")
    en = 0.0
    for i in range(len(k)):
        en = en - epsilon(k[i],h,j)
    en = en/numpy.float(L)
    return en

# coefficients of the products that gives c_\alpha
def Coef( h0, h, L ):
    j = 1.0
    k = kpoints(L,0,">")
    c = []
    for i in range(len(k)):
        ek  = epsilon( k[i], h, j )
        e0k = epsilon( k[i], h0, j )
        c.append( (4.0*(h - h0)*(h - h0) - (e0k - ek)*(e0k - ek))/(-4.0*(h - h0)*(h - h0) + (e0k + ek)*(e0k + ek)) )

    return numpy.array(c)

def CoefGS( h0, h, L ):
    j = 1.0
    k = kpoints(L,0,">")
    proj = 1.0
    ek  = epsilon( k, h, j )
    e0k = epsilon( k, h0, j )

    for i in range(len(k)):
        proj = proj*(-4.0*(h - h0)*(h - h0) + (e0k[i] + ek[i])*(e0k[i] + ek[i]) )/( 4.0 * ek[i] * e0k[i] )

    return proj
