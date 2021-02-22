'''Laxman Pandey, Fall 2015'''
from __future__ import division
import matplotlib.pylab as plt
import numpy as np

pi = np.pi
sin = np.sin
cos = np.cos
sqrt = np.sqrt
h = 6.626e-34 # Js
hbar = h/(2*pi)
me = 9.10938291e-31 # kg
c = 299792458 # m/s
fact = np.math.factorial

# Particle in a box
def pib_wavefunc(n,L):
    x = np.linspace(0,L,500)
    y_n = sqrt(2/L)*sin(n*pi*x/L)
    return y_n, x

def pib_plot(n2,L,n1=1,wavefuncPlot=True,probDenPlot=False):# wavefunction plots
    E_nim1 = 0
    for ni in range(n1,n2+1):    
        E_ni = ni**2*h**2/(8*me*L**2)
        deltaE = E_ni - E_nim1
        print('E_%s = %.2e J  deltaE = %.2e J' %(ni, E_ni, deltaE))
        y_ni, xi = pib_wavefunc(ni, L)
        if wavefuncPlot:
            plt.plot(xi, y_ni, label=r'$\psi_%s$' %ni)
            #plt.ylabel(r'$\psi_%s$' %ni)
            plt.ylabel('Amplitude')
            E_nim1 = E_ni

        if probDenPlot:
            plt.plot(xi, y_ni**2, label=r'$\psi_%s^*\psi_%s$' %(ni,ni))
            plt.ylabel('Amplitude')
        plt.xlabel('L')
        plt.title('Particle In a Box')
        plt.legend()
    plt.show()

# Particle in a ring
        
def pir_wavefunc(n,phi):
    ang = np.linspace(0,phi,500)
    y_n = (1/sqrt(2*pi))*np.exp(1j*n*phi)
    return y_n, ang

def pirE(n,R):
    return n**2*h**2/(8*pi**2*me*R**2)

def pirDeltaE(n1,n2,R):
    return pirE(n2,R) - pirE(n1,R)

# For benzene using PIR
print('Transition Energy (n1-n2 PIR)%.2e' %pirDeltaE(1,2,1.39e-12))
print('%.2e' %(h*c/pirDeltaE(1,2,1.39e-12)))

# Harmonic Oscillator
def hermite(x,n):
    '''returns nth hermite polynomial'''
    if n==0:
        return 1
    elif n==1:
        return 2*x
    else:
        return 2*x*hermite(x,n-1) - 2*(n-1)*hermite(x,n-2)

def ho_wavefunc(n,a):
    '''returns the simple harmonic oscillator wavefunction for n = 0, 1, 2, ...
    a = alpha*x             # a is dimensionless 
    alpha = (m*omega/hbar)**0.5
    omega = (k/m)**0.5      # angluar frequency
    m = reduced mass
    '''
    y_n = (1/(pi**0.5*2**n*fact(n)))**0.5*hermite(a,n)*np.exp(-a**2/2)
    return y_n, a

def ho_plot(n2,alpha,n1=0,wavefuncPlot=True,probDenPlot=False):
    '''wavefunction and probability density plots
    n1 = smallest quantum number
    n2 = largest quntum number
    '''
    x0 = ((2*n2+1)*hbar/(m*omega))**0.5
    x = np.linspace(-2*x0,2*x0,1000)
    a = alpha*x
    for ni in range(n1,n2+1):
        y_ni, xi = ho_wavefunc(ni,a)
        if wavefuncPlot:
            plt.plot(xi, y_ni, label=r'$\psi_%s$' %ni)
            #plt.ylabel(r'$\psi_%s$' %ni)
            plt.ylabel('Amplitude')

        if probDenPlot:
            plt.plot(xi, y_ni**2, '--', label=r'$\psi_%s^*\psi_%s$' %(ni,ni))
            plt.ylabel('Amplitude')
        plt.xlabel('x')
        plt.title('Harmonic Oscillator')
        plt.legend()
    plt.show()
    #plt.savefig('ho_psi0.pdf')

# call functions
pib_plot(2,4,2,wavefuncPlot=True,probDenPlot=True) # plot wavefunctions, prob. densities, print energies
m = 1.6267E-27 #kg
omega = 5.63212E14 # rad/s
alpha = (m*omega/hbar)**0.5
ho_plot(0,alpha,0,wavefuncPlot=True,probDenPlot=True)


# Hayward Problem 2.4
L = 22*144e-12 # m
a = (144/2 + 10*144)*1e-12 # m # starting halfway 10-11
b = a + 144e-12 # m # 1 bond after starting point
n = range(1,12) # 1, 2, 3, ..., 11
ni, nf = 11, 12

def deltaE(ni,nf):
    dE = (nf**2 - ni**2)*h**2/(8*me*L**2)
    return dE

dE = deltaE(ni,nf)
def print_deltaE(dE):
    print('deltaE = %.2e J' % dE)
    print('absorbance wavelength = %.4e m' %(h*c/dE))

print_deltaE(dE)


# Calculating probability of finding a particle between x=a and x=b
'''
P_a2b = (2/L) * Integral of sin^2(n*pi*x/L) dx from a to b.
      = (1/L)*(b-a) - 1/(2*n*pi) * sin(2*n*pi*b/L) + 1/(2*n*pi) * sin(2*n*pi*a/L)
'''

def pib_prob_a2b(n,L,a,b):
    '''Returns probability for finding a particle between x=a and x=b for a PIB model
    n: quantum number, an integer >= 1
    L: length of the box
    a and b: lower and upper limits of the range/boundary'''
    P_a2b = (1/L)*(b-a) - 1/(2*n*pi) * sin(2*n*pi*b/L) + 1/(2*n*pi) * sin(2*n*pi*a/L)
    return P_a2b

# Using Approximation as suggested in the Hayward problem
'''
P_a2b = average of (2/L) * ( (sin^2(n*pi*a/L) + sin^2(n*pi*b/L) )
'''
width = 144e-12 #m width of the range where approximated
def pib_prob_a2b_approx(width,n,L,a,b):
    P_a2b_approx = width*0.5*((2/L)*(sin(n*pi*a/L))**2 + (2/L)*(sin(n*pi*b/L))**2)
    return P_a2b_approx

# compute total probabilities
def print_tot_probability(exact=True,approx=True):
    P_a2b, P_a2b_approx = 0, 0 # set null probabilities
    for ni in n:
        P_a2b += pib_prob_a2b(ni,L,a,b)
        P_a2b_approx += pib_prob_a2b_approx(width,ni,L,a,b)

    if exact:
        print('\nAvg. # of pi el. betwn. %.2e to %.2e for L=%.3e: %.2f (exact)' % (a,b,L,P_a2b * 2)) # *2 for 2 e/level
    if approx:
        print('Avg. # of pi el. betwn. %.2e to %.2e for L=%.3e: %.2f (approx.)' % (a,b,L,P_a2b_approx * 2)) # *2 for 2 e/level
    pass
print_tot_probability(exact=True,approx=True) # print total probabilities


# Homework 2 Problems
L = 12*144e-12 # m; box of length 12*144 pm
#a, b, width = (1/3)*L, (2/3)*L, (1/3)*L # central third of the box
a, b, width = (5/12)*L, (6/12)*L, (1/12)*L # small width of the box
n = range(1,6) # n = 1, 2, 3, ..., 6

def print_probability(exact=True,approx=False):
    for ni in n: # for 
        if exact:
            print('\nHW2.2-3: For n=%d, prob. of e from %.2e to %.2e for L=%.3e: %.3f (exact)' % (ni,a,b,L,pib_prob_a2b(ni,L,a,b)))
        if approx:
            print('HW2.4:   For n=%d, prob. of e from %.2e to %.2e for L=%.3e: %.3f (approx.)' % (ni,a,b,L,pib_prob_a2b_approx(width,ni,L,a,b)))
        pass
print_probability(exact=True,approx=True) # print probabilities

dE = deltaE(6,7) # HOMO to LUMO
print_deltaE(dE)
print_tot_probability(exact=True,approx=True) # print total probabilities



