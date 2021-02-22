#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      lpandey
#
# Created:     05/09/2014
# Copyright:   (c) lpandey 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------
from __future__ import division
#from sys import exit
import math as mt
import re
import numpy as np
from matplotlib import pyplot as plt

# elements
atomicMasses = {
'H': 1.00797, 'He': 4.0026, 'Li': 6.939, 'Be': 9.0122, 'B': 10.811,
'C': 12.0112, 'N': 14.0067, 'O': 15.9994, 'F': 18.9984, 'Ne': 20.183, 
'Na': 22.9898, 'Mg': 24.312, 'Al': 26.9815, 'Si': 28.086, 'P': 30.9738, 
'S': 32.064, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.102, 'Ca': 40.08, 
'Sc': 44.956, 'Ti': 47.9, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938, 
'Fe': 55.847, 'Co': 58.9332, 'Ni': 58.71, 'Cu': 63.54, 'Zn': 65.37, 
'Ga': 69.72, 'Ge': 72.59, 'As': 74.9216, 'Se': 78.96, 'Br': 79.909, 
'Kr': 83.8, 'Ba': 137.34, 'I': 126.904, 'Pb': 207.19, 'Cs': 132.905, 
'Ag': 107.87, 'Sn': 118.69, 'Au': 196.967, 'Hg': 200.59, 'Xe': 131.3        
}

def getAtomicMass(el):
    return atomicMasses['%s' % el]

def getFormulaMass(formula):
    '''
    Returns total mass of elements in a formula, eg. 'Ca(NO3)2'
    '''
    formulaMass = 0
    atomslist = atoms(formula)
    #print(type(atomslist), len(atomslist), atomslist)
    for k, v in atomslist.pop().items():
        formulaMass += (getAtomicMass(k) * int(v))
    print('%s: %.2f g/mol' %(formula, formulaMass))
    return formulaMass

def atoms(formula,debug=False,stack=[], delim=0,
    atom = r'([A-Z][a-z]?)(\d+)?', ldel=r'\(', rdel=r'\)(\d+)?'):
    '''
    Returns formula mass of the molecule Eg. 'CH3COOH' 
    '''
    formula = formula.replace(' ', '')
    formula = str(formula)
    #atom = r'([A-Z][a-z]?)(\d+)?' #regex for string matching atom
    #ldel=r'\(' #regex for '('
    #rdel=r'\)(\d+)?' #regex for ')'
    #stack = []
    
    re_atom = re.match(atom, formula) #match object
    #print(formula, re_atom)
    re_ldel = re.match(ldel, formula)
    re_rdel = re.match(rdel, formula)
  
    if re_atom:
        tail = formula[len(re_atom.group()):] #remaining atoms
        head = re_atom.group(1) #matched atom
        num  = int(re_atom.group(2) or 1) #repetition if exists else 1 
        stack = stack or [{}]      
        stack[-1][head] = stack[-1].get(head, 0) + num #value (repetition) for key (atom)
        if debug:
            print([head, num, tail])
    
    elif re_ldel:
        tail   = formula[len(re_ldel.group()):] #minus '('
        delim += 1
        stack.append({})
        if debug:
            print(['left-delimiter', tail])
    
    elif re_rdel:
        tail   = formula[len(re_rdel.group()):] 
        num    = int(re_rdel.group(1) or 1) #repitition number
        delim -= 1
        if delim < 0:
            raise SyntaxError("un-matched right parenthesis in '%s'" % (formula,))
        #for (k, v) in stack.pop().iteritems():#python 2.x
        for (k, v) in stack.pop().items():#python 3.x
            stack[-1][k] = stack[-1].get(k, 0) + v*num
        if debug:
            print(['right-delimiter', num, tail])

    else:
        raise SyntaxError("'%s' does not match any regex" % (formula,))

    if len(tail) > 0:
        #print(tail)
        atoms(tail, debug, stack, delim, atom, ldel, rdel)
        return stack      
    else:
        if delim > 0:
            raise SyntaxError("un-matched left parenthesis in '%s'" % (formula,))
        if debug:
            print(stack[-1])
        return stack

#Physical
# define constants
c =299792458 # m/s
h = 6.62606957e-34 # m2 kg / s
kB =1.3806488e-23 # m2 kg s-2 K-1 # J/K
kBw = 0.6950 # cm-1
NA =6.0221413e+23 # particles
amu = 1.66053892e-27 # kg
R = 8.314 # J /molK
RLatm = 8.206e-2 #L atm/molK
g = 9.80665 # m/s2
permitivity = 8.85418782e-12 # C V-1 m-1
me = 9.10938291e-31 # kg
mp = 1.6726219e-27 # kg
e = 1.60217657e-19 # coulombs
RH = (me*e**4/(8*permitivity**2*h**2)) # J
pi = mt.pi
def ln(n):
    return np.log(n)
def hartrees2eV(Eh):
    return Eh*27.2114 # eV
def hartrees2kcalmol(Eh):
    return Eh*627.503 # kcal/mol
def invCm2eV(invCm):
  return float(invCm)*1.23984e-4 # eV
def invCm2J(invCm):
  return float(h*c*invCm*1e2) # J
def invCm2K(invCm):
    return invCm2eV(invCm)/8.621738e-5
hbar = h/2/pi
br = 0.52917721067e-10 # bohr radius in m
def exp(n):
    return np.exp(n)
def invCm2hartrees(invCm):
  return float(invCm)*4.55633e-6 # Eh
def J2eV(J):
  return float(J)*6.24150934e18 # eV
def J2nm(J):
  return h*c*1e9/float(J) # m
def eV2J(eV):
  return eV*e # J
def K2eV(K):
  return J2eV((K)*kB) # eV
def K2invCm(K):
  return K2eV(K)*8065.5 # cm-1
def nm2eV(nm):
  return J2eV(h*c/(float(nm)*1e-9))
def C2K(C):
  return float(C)+273.15 # K
def K2C(K):
  return float(K)-273.15 # C
def deg2rad(dg):
    return 2*mt.pi*dg/360
def rad2deg(rad):
    return rad*180/mt.pi
def wavenum2nm(invCm):
    return (1/invCm)*1e7
def hydrogenWavenum(n1, n2=float('nan')):
    if mt.isnan(n2):
        return RH*(1/n1**2)
    else:
        return RH*(1/n1**2 - 1/n2**2)
def BohrRadius(n,Z):
   """returns Bohr radius in m given n and Z"""
   return (4*pi*permitivity*n**2*hbar**2)/(me*e**2*Z) # in m
def BohrESpeed(n,Z):
   """returns electron speed in Bohr orbit (m/s) given n and Z"""
   return (e**2*Z)/(4*pi*permitivity*n*hbar) # in m
def BohrEnergy(n,Z):
    """returns Bohr Energy in J given n and Z"""
    return (-me*e**4*Z**2/(8*permitivity**2*h**2*n**2)) # in eV
def RydbergTransition(nf,ni,Z):
    """returns photon wavelength (nm) corresponding to transition from
    ni to nf in Bohr atom for given Z"""
    return abs(J2nm(BohrEnergy(nf,Z)-BohrEnergy(ni,Z)))

def rdms(m1,m2): # reduced mass
    return m1*m2/(m1+m2)

def amu2kg(m): # amu to kg
    return m*amu #

def pH(hydronium):
    """returns  negative log of hydronium ion concentration
    """
    return -np.log10(hydronium)

def MaxwellDistribution(M,T,v):
    """returns the Maxwell distribution function for M, T, v"""
    return 4*pi*(M/(2*pi*R*T))**(3/2)*(v**2)*exp(-M*v**2/(2*R*T))

def speed_rms(M,T):
    """returns rms speed of particles with mass M and temp T"""
    return (3*R*T/M)**0.5

def speed_mean(M,T):
    """returns mean speed of particles with mass M and temp T"""
    return (8*R*T/(pi*M))**0.5

def speed_mp(M,T):
    """returns most probable speed of particles with mass M and temp T"""
    return (2*R*T/M)**0.5

def pibE(n=1,m=mp,L=1e-10,d=1):
    if d==3:
        nx, ny, nz = [int(x) for x in input('Enter states nx ny nz: ').split()]
        a, b, c = [float(x) for x in input('Enter sides a b c: ').split()]
        return (nx*h)**2/(8*m*a**2)*(ny*h)**2/(8*m*b**2)*(nz*h)**2/(8*m*c**2)
    elif d==1:
        return (n*h)**2/(8*m*L**2)
    
def flowRate(r,n,dP,l):
    """returns Poiseuille's flow rate given r, n, dP and l"""
    return (pi/8)*(r**4/n)*(dP)*(1/l)

def FlowRatio(r1,r2,dP1,dP2):
    return (r1/r2)**4 * (dP1/dP2)

def percentError(expt,theor):
    '''returns percent error when experimental(expt) and
    theoretical(theor) values are provided'''
    return 100*((expt-theor)/theor)

def quadratic(a,b,c):
    """a, b, and c are coefficients of x^2, x^1, and x^0.
    For example, 3x^2 + 4x -2 would have ...
    ... a = 3, b = 4, and c = -2.
    """
    try:
        x1 = (-b + mt.sqrt(b**2 - 4*a*c))/(2*a)
        x2 = (-b - mt.sqrt(b**2 - 4*a*c))/(2*a)
        return x1, x2
    except ValueError:
        return 'At least one root is imaginary'

def coord2distance(pt1, pt2=0):
    #points are entered as coordinates (x,y,z)
    if pt2 == 0:
        dx, dy, dz = pt1[0], pt1[1], pt1[2]
    else:
        dx = pt1[0]-pt2[0]
        dy = pt1[1]-pt2[1]
        dz = pt1[2]-pt2[2]
    return mt.sqrt(dx**2 + dy**2 + dz**2)

def coord2angle(pt1, pt2, pt3=0, pt4=0):
    #points are entered as coordinates (x,y,z)
    # vectors made from pt2 to edge points pt1 and pt3
    v12 = np.array(pt2) - np.array(pt1)
    if pt3 == 0:
        v13 = -np.array(pt1)
        v23 = -np.array(pt2)
    else:
        v13 = np.array(pt3) -np.array(pt1)
        v23 = np.array(pt3) -np.array(pt2)
    # normalized v1 and v3
    # vector/length of vector
    n12 = v12/np.linalg.norm(v12)
    n13 = v13/np.linalg.norm(v13)
    n23 = v23/np.linalg.norm(v23)

    #angles
    a123 = np.arccos(np.dot(n12,n23))
    a132 = np.arccos(np.dot(n13,n23))
    a312 = np.arccos(np.dot(n13,n12))
    return rad2deg(a312), rad2deg(a123), rad2deg(a132)

def coord2dihedral(pt1, pt2, pt3, pt4=0):
    #points are entered as coordinates (x,y,z)
    d12 = coord2distance(pt1, pt2)
    d23 = coord2distance(pt2, pt3)
    d34 = coord2distance(pt3, pt4)
    d14 = coord2distance(pt1, pt4)
    d24 = coord2distance(pt2, pt4)
    d13 = coord2distance(pt1, pt3)

    p = d12**2 * ( d23**2 + d34**2 - d24**2) + \
        d23**2 * (-d23**2 + d34**2 + d24**2) + \
        d13**2 * ( d23**2 - d34**2 + d24**2) - \
        2 * d23**2 * d14**2

    q = (d12 + d23 + d13) * ( d12 + d23 - d13) * \
        (d12 - d23 + d13) * (-d12 + d23 + d13) * \
        (d23 + d34 + d24) * ( d23 + d34 - d24) * \
        (d23 - d34 + d24) * (-d23 + d34 + d24)

    da = np.arccos(p/np.sqrt(q))
    return rad2deg(da)


def deBroglie()
    return h/p

def permutation(n,j):
    '''returns number of ordered subset of j objects from a total of n'''
    return mt.factorial(n)/mt.factorial(n-j)

def combination(n,j):
    '''returns number of unordered subset of j objects from a total of n'''
    return mt.factorial(n)/mt.factorial(n-j)/mt.factorial(j)

#plot Boyle's Law
def BoylesLaw(constant,lowP,highP,plotType='PV'):
    '''plots Boyles Law
    provided constant and low and high pressure or volumes'''
    ps = np.linspace(lowP,highP,200)
    if plotType == 'PV' or plotType == 'pv':
        vs = constant/ps
        plt.plot(vs,ps)
        plt.xlabel('Volume ($m^3$)')
    
    else:
        invV = ps/constant
        plt.plot(invV,ps)
        plt.xlabel('1 / Volume ($m^{-3}$)')
    plt.ylabel('Pressure (Pa)')
    plt.show()

#BoylesLaw(10,1,10,plotType='PiV')
#BoylesLaw(10,1,10,plotType='PV')

##print('%e' %(13.6*100**3/1000))
##print('%s' %(750*2*getTotalMass('Fe')/getTotalMass('Fe2 O3')))
##print('%e J' %(h*c/564e-9))
##print('%e eV' %(J2eV(h*c/564e-9)))
##print('%.3g' % (1.23*(1/getTotalMass('N O O'))))
##print('%.2e' % ((0.917/1)*(1/1000)*(1/1000)*(2.54**3/1)*(12**3/1)*(5.5e6)*(5280**2/1)*(6000/1)))
#print 'Q3: %.3e' %(1.8e3*1.01e5/14.7)
#print 'Q5: %.3e' %(2.30*0.08206*C2K(23)/(getTotalMass('C3 H8')*0.250))
#print 'Q7: %.3e' %(2.49*0.08206*C2K(62)/(1.98*0.752))
#print 'Q10: %.3e' %(4.00*50*10.0/(0.08206*C2K(19)))
#print 'Q15: %.3e' %(1/(4e-4*600+1/0.8))
#print 'Q17: %.3e' %(0.991*0.08206*298.15/1.2)
#print 'Q21: %.3e' %(100**2/(500*30**3))
#print 'Q22: %.3e' %(0.55**2/(0.25*0.035))
#print 'Q100: %.3e' %(NA*25/22.4e3)
#print 'Q101: %.3e' %(22.4*4.5e18/NA)
#print 'Q102: %.3e' %(4.5e18/NA)
#print 'Q103: %.3e' %(9144/115)
#print 'Q104: %.3e' %(3e10*60*60/(2.54*12*5280))
#print 'Q105: %.3e' %(0.075e3/58.5)
#print 'Q106: %.3e' %(1.5e-3*58.5e2)
#print('Q106: %.3e' %(2.33e-4*3.25e10/7))

def main():
    pass

if __name__ == '__main__':
    main()
