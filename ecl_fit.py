#!/usr/bin/python
#Filename: ecl_fit.py

# Evaluate Effective Conjugation Length for pi-conjugated oligomeric
# systems using several fit methods ...
# Meier (exponential) fit (Acta Polymer 1997)
# Modified Kuhn Fit (Gierschner, Adv. Mater. 2007)
# Modified Kuhn Fit (Slava, Research Scientist in our Group)

from __future__ import division

import sys
# get data files
if len(sys.argv) < 2:
  print "\n Usage: ecl_fit.py datafiles\n (eg. ecl_fit.py file1.data file2.data)\n"
  sys.exit(1)

import os
import numpy
from scipy.optimize import leastsq
from scipy import polyfit, stats
from pylab import *
from matplotlib.font_manager import FontProperties

# packages needed : dependencies, numpy, scipy, pylab (or matplotlib)

EPSILON = 1e-8 # Smallest significant difference
MAXFEV = 50000 # Max function evaluation

# data files
dfs = sys.argv[1:]
#print '\nNumber of samples : %s \n' % len(dfs)

# get data function
def getData(f):
  if not (os.path.isfile('%s' % f)):
    print ' %s DOES NOT EXIST !!!' % f
    sys.exit(1)
  n, nm = [], [] # assuming data in nanometers
  rf = open('%s' %f, 'r')
  lines = rf.readlines()
  for line in lines:
    line = line.strip()
    if not line: continue
    # checking Units
    if ' nm' in line: unit = 'nm'
    elif ' cm-1' in line: unit = 'cm-1'
    elif ' eV' in line: unit = 'eV'
    line = line.split()
    if len(line) < 2: continue
    if line[0][0] == '#': continue # ignore comment line
    try:
      ln, lnm = eval(line[0]), eval(line[1])
      if type(ln) == int: n.append(ln)
      if (type(lnm) == int or type(lnm) == float): nm.append(lnm)
    except: continue
  rf.close()
  print "\n%s read units (%s) write units (%s)" % (f, unit, 'eV')
  # Unit conversions to eV  
  if unit == 'eV':
    return array(n), array(nm)     # already in eV 
  elif unit == 'nm':
    return array(n), 1240/array(nm) # convert from nm to eV
  elif unit == 'cm-1':
    #return array(n), array(nm)*0.00012398 # convert from cm-1 to eV
    return array(n), array(nm)*0.00012398*1000 # convert from x1000 cm-1 to eV

# now call function to get data
for df in dfs:
  dn, dEn = getData(df)
  assert len(dn) == len(dEn), '%s %s %s' % (df, len(dn), len(dEn))
  #plot(dn, dEn, 'o', markersize=10)
  plot(1/dn, dEn, 'o', markersize=10)

  # WHEN EXCLUDING MONOMER DATA
  dn, dEn = dn[1:], dEn[1:]
  E1 = dEn[0]

  # polynomial fit
  # n[i:] used to include or exclude monomer data
  def select(i, degree):
    linearls = []
    dn1 = dn[i:]; dEn1 = dEn[i:]

    # polynomial regression
    coeffs = numpy.polyfit(1/dn1, dEn1, degree)
    # polynomial coefficients
    results = coeffs.tolist()
    # r-squared
    p = numpy.poly1d(coeffs)
    yhat = [p(z) for z in 1/dn1]
    ybar = sum(dEn1)/len(dEn1)
    ssreg = sum([(yihat - ybar)**2 for yihat in yhat])
    sstot = sum([(yi - ybar)**2 for yi in dEn1])
    rsq = ssreg / sstot
    # this is assuming linear function, degree 1
    nx = linspace((min(dn1)), 110, num=2000)
    #print dn1, dEn1
    m1,b1 = polyfit(1/dn1, dEn1, degree)
    plot(1/nx, b1 + m1/nx, '--', linewidth=3, label='Linear')
    ny = linspace((min(dn1)), 100, num=400)
    print "%3s %7s %7s %9s" % ('n', 'eV', 'nm', 'cm-1')
    for i, j in zip(dn1, dEn1):
      print '%3d %7.3f %7.1f %9.1f' % (i, j, 1240/j, j*8065.5)
    for ik in range(len(ny)):
      linearls.append('%9.5f %9.5f' %(ny[ik], b1 + m1/ny[ik]))
    return results[1], rsq, linearls

  # Now call function select &
  # collect linear (pylynomial) fit data
  lin_Einf, lin_rsq, lin_ls = select(0,1) 
  #lin_Einf, lin_rsq, lin_ls = select(1,1) 

  # other fits follow
  # define Meier (exponential) fit (Acta Polymer 1997)
  def ExpnEn(n, Einf, a):
    return Einf + (E1 - Einf) * exp(-a*(n-1))

  # define Modified Kuhn Fit (Gierschner, Adv. Mater. 2007)
  def KuhnEn(n, A, B):
    return numpy.sqrt(A + B * numpy.cos(numpy.pi / (n + 1.0))) 

  def KuhnEn2(n, B):
    return numpy.sqrt(E1**2 + B * numpy.cos(numpy.pi / (n + 1.0)))

  # define Modified Kuhn Fit (Slava)
  def SlavaEn(n, C, D):
    return C + D * numpy.cos(numpy.pi / (n + 1.0))

  def SlavaEn2(n, D):
    return E1 + D * numpy.cos(numpy.pi / (n + 1.0))

  # initial values for paramters
  init = 1.0       # initial parameter value

  # intial undetermined parameters
  ks0 = [init, init]

  # called with an array of undetermined parameters

  def residualsE(ks):
    Einf, a = ks
    y = ExpnEn(dn[n_i::], Einf, a)
    return y - dEn[n_i::]

  def residualsK(ks):
    A, B = ks
    y = KuhnEn(dn[n_i::], A, B)
    return y - dEn[n_i::]

  def residualsK2(ks):
    B = ks
    y = KuhnEn2(dn[n_i::], B)
    return y - dEn[n_i::]

  def residualsS(ks):
    A, B = ks
    y = SlavaEn(dn[n_i::], A, B)
    return y - dEn[n_i::]

  def residualsS2(ks):
    B = ks
    y = SlavaEn2(dn[n_i::], B)
    return y - dEn[n_i::]

  Expnls, Kuhnls, Slavals = [], [], []

  for n_i in [0]:
  #for n_i in [0,1]: # can have 0 excluded by slicing out 0    
    def getFit(myresidual, ks0T):
      ksT,_ = leastsq(myresidual, ks0T,
                      maxfev=MAXFEV, ftol=EPSILON, xtol=EPSILON)
      return ksT

    # Now call funtion myks to get the fitted parameters
    Einf, a,  = getFit(residualsE, ks0)
    A, B = getFit(residualsK, ks0)
    B2 = getFit(residualsK2, init)
    C, D = getFit(residualsS, ks0)
    D2 = getFit(residualsS2, init)

    Kuhn_Einf = numpy.sqrt(A + B)
    Kuhn_Einf2 = numpy.sqrt(E1**2 + B2)
    Slava_Einf = C + D   
    Slava_Einf2 = E1 + D2   

    nx = linspace((min(dn)), 110, num=2000)
    ny = linspace((min(dn)), 100, num=400)

    #print dn, dEn
    for ik in range(len(ny)):
      Expnls.append('%9.5f %9.5f' %(ny[ik], ExpnEn(ny[ik], Einf, a)))
      Kuhnls.append('%9.5f %9.5f' %(ny[ik], KuhnEn(ny[ik], A, B)))
    #plot(1/nx, ExpnEn(nx, Einf, a), '-', linewidth=3, label='%s' % 'Exp.')
    plot(1/nx, KuhnEn(nx, A, B), '-', linewidth=3, label='%s' % 'Kuhn')
    plot(1/nx, KuhnEn2(nx, B2), '-', linewidth=3, label='%s' % 'Kuhn2')
    plot(1/nx, SlavaEn(nx, C, D), '-', linewidth=3, label='%s' % 'Slava')
    plot(1/nx, SlavaEn2(nx, D2), '-', linewidth=3, label='%s' % 'Slava2')

    #print "\nFit Parameters : Einf R**2 a A and C where appropriate"
    print "\nFit Parameters : Einf R**2 A and C where appropriate"
    print "Linear  : %.3f %.3f" % (lin_Einf, lin_rsq)
    print "Expn    : %.3f %.3f" % (Einf, a)
    print "Kuhn    : %.3f %.3f" % (Kuhn_Einf, numpy.sqrt(A))
    print "Kuhn2   : %.3f"      % (Kuhn_Einf2)
    print "Slava   : %.3f %.3f" % (Slava_Einf, C)
    print "Slava2  : %.3f"      % (Slava_Einf2)
            
  # plot editing
  #title('First excitation energy as a function of inverse repeat units', 
  #       fontsize=16)
  l = legend(loc='upper left', prop={'size':20, 'weight':'bold'})
  #l = legend(loc='upper center')
  l.draw_frame(False)
  #ylim(1.00, 3.50)
  xlim(0.00, 1.05)
  #xlim(0.00, 0.55)
  if df == 'Nakanishi1998_Oct.data':
    xlim(0.00, 0.2)

  xl = xlabel('1/n', fontsize=20, fontweight='bold')
  yl = ylabel('First transition energy (eV)', fontsize=20, fontweight='bold')
  #setp(gca(), 'yticklabels', 
  
  #write data for gnulpot
  wf = open("%s.fit" % df.replace('.data', ''), 'w')
  wf.write('%20s%20s%20s\n' 
            %('#linearFit (n, E1)', '#ExpnFit (n, E1)', '#KuhnFit (n, E1)'))
  assert len(Expnls) == len(Kuhnls) == len(lin_ls)
  for dt in range(len(Expnls)):
    wf.write('%s   %s   %s\n' % (lin_ls[dt], Expnls[dt], Kuhnls[dt])) 
  wf.close()

  savefig('%s.png' % df.replace('.data', ''))
  print "###########################"
show()
    
#from time import sleep
#sleep(5)
