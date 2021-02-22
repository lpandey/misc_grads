#!/usr/bin/env python

import sys
from numpy import sqrt, arange, log, exp, pi

if len(sys.argv) < 2:
  print '\nUsage: gaussian_broadening.py inputfile'
  print ' (eg. gaussian_broadening.py file.dat44\n'
  sys.exit(1)

ifile=sys.argv[1]
try:
  rf = open(ifile, 'r')
except:
  print 'Cannot open %s' % ifile
  sys.exit(1)

NT = int(rf.readline().replace('#', ''))
line2 = rf.readline().replace('#', '').split()
FWHM, Emin, Emax = float(line2[0]), float(line2[1]), float(line2[2])

count=0
Es, fs = [], []
lines=rf.readlines()
rf.close()
for line in lines:
  line = line.strip()
  if not line: continue
  line = line.split()
  E, f = float(line[0]), float(line[1])
  Es.append(E)
  fs.append(f)
  count += 1

assert NT == count, 'Number of states mismatch'
wf = open('%s55' % ifile[:-1], 'w')

dE = (Emax - Emin)/3000.0
FWHM2 = FWHM**2/4.0
for i in arange(3001):
  Ex    = Emin + (i)*dE
  Wlx   = 1240.00/Ex
  invCm = 8065.73*Ex
  Yx    = 0.0
  for j in arange(NT):
    ER  = Ex - Es[j]
    ERG = log(2.0)*ER**2/FWHM2 
    Exx = fs[j]/exp(ERG)
    Yxx = Yx + Exx
    TMP = 0.5*sqrt(pi/log(2))  
    Abs = Yxx/(4.32E-9*FWHM*8065.73*0.5*TMP)
  #print Ex, invCm, Wlx, Yxx, Abs
  wf.write('%9.3f%10.1f%10.2f%7.2f%13.1f\n' % (Ex, invCm, Wlx, Yxx, Abs))

wf.close()

