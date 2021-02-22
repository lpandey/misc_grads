#! /usr/bin/env python
# Laxman Pandey 5.3.2012

# Calculate Coulomb Energy of two molecules based on charges on atoms
# for each molecule and their separation.
# E(Coul)=Sumover[atomD, atomA] {q(atomA) * q(atomB) / permitivity*r(atomA,atomD)}  
import os
import sys
from math import sqrt, pi
from subprocess import Popen, PIPE
from itertools import islice

if len(sys.argv) < 3:
    print "\n Usage: coulombEnergy_J.py NwchemOutFile.out numAtomsFrag2"
    print " (eg. coulombEnergy_J.py p3ht_c1_m1.out 98)\n" 
    sys.exit(1)

logf=sys.argv[1]; nAtoms2=sys.argv[2]
try:
  nAtoms2=eval(nAtoms2)
except:
  print "\nPlease review number of atoms for the second Fragment\n"   
  sys.exit(1)

if not os.path.exists(logf):
  print "\n %s Does NOT exist.\n" % logf
  sys.exit(1)

# define grep commands
geomLnCmd='grep -nm 1 "%s" %s' %('XYZ format geometry',logf)
chkCmd='grep -c "%s" %s' %('Total times  cpu', logf)
eLnCmd='grep -nm 2 "%s" %s' %('Total Density - Mulliken Population Analysis', logf)

# Check normal job termination
chk=geomLn=Popen('%s' %chkCmd, shell=True, stdout=PIPE).stdout
chk=chk.read().strip()
if int(chk)==0:
  print ' !! %s did not terminate normally !!' %logf 
  #sys.exit(1)

# define functions
def seek_to_line(f, n):
  for ignored_lines in islice(f, n):
    pass # skip beginning n lines

def myCoords(myf):
  fragA,fragB=[],[]
  geomLn=Popen('%s' %geomLnCmd, shell=True, stdout=PIPE).stdout
  geomLn=geomLn.read().strip().split()[0][:-1]
  #print geomLn
  rf=open(myf, 'r')
  seek_to_line(rf, int(geomLn))
  lcount=0
  for line in rf:
    line=line.strip()
    if not line: continue
    line=line.split()
    if len(line)==1:
      try:
        nAtomsTot=eval(line[0])
      except: continue
    nAtoms1=nAtomsTot-nAtoms2
    if len(line)!=4: continue 
    el,x,y,z=line[0],line[1],line[2],line[3]
    assert len(el)<=2, 'is %s an element?' % el
    assert type(eval(x))==float, 'is %s x-coordinate?' % x
    assert type(eval(y))==float, 'is %s y-coordinate?' % y
    assert type(eval(z))==float, 'is %s z-coordinate?' % z
    lcount+=1
    if lcount <= nAtoms1:
      fragA.append('%s %s %s %s' %(el,x,y,z))
    else:
      fragB.append('%s %s %s %s' %(el,x,y,z))
    #print line, lcount
    if lcount==nAtomsTot: break
  rf.close()
  return nAtomsTot, fragA, fragB

nAtomsTot,cA,cB=myCoords(logf)
nAtoms1=nAtomsTot-nAtoms2
#print cA
#print cB

def myCharges(myf):
  efragA,efragB=[],[]
  eLn=Popen('%s' %eLnCmd, shell=True, stdout=PIPE).stdout
  eLn=eLn.readlines()[-1:][0].split()[0][:-1]
  #print eLn
  rf=open(myf, 'r')
  seek_to_line(rf, int(eLn))
  lcount=0
  for line in rf:
    line=line.strip()
    if not line: continue
    line=line.split()
    if len(line) < 7: continue
    n,el,c0,c=line[0],line[1],line[2],line[3]
    assert type(eval(n))==int, 'is %s atom number?' % n
    assert len(el)<=2, 'is %s an element?' % el
    assert type(eval(c0))==int, 'is %s atomic number?' % c0
    assert type(eval(c))==float, 'is %s charge?' % c
    lcount+=1
    if lcount <= nAtoms1:
      #efragA.append('%s %s' %(el,float(c))) # g09 format
      efragA.append('%s %s' %(el,float(c)-int(c0))) # nwchem format
    else:
      #efragB.append('%s %s' %(el,float(c))) # g09 format
      efragB.append('%s %s' %(el,float(c)-int(c0))) # nwchem format
    #print line, lcount
    if lcount>=nAtomsTot: break
  assert(int(n)==nAtomsTot), '!! ALERT atom number check FAILED !!'
  rf.close()
  return efragA, efragB

eA,eB=myCharges(logf)
#print eA
#print eB

# Define function to calculate Coulomb energy
def myCoulomb(fragA,fragB,efragA,efragB):
  myJeV,myJJ=0,0
  for i in range(len(fragA)):
    for j in range(len(fragB)):
      coordsA=fragA[i].split()
      xA,yA,zA=coordsA[1],coordsA[2],coordsA[3]
      xA,yA,zA=float(xA),float(yA),float(zA)
      qA=float(efragA[i].split()[1])
      coordsB=fragB[j].split()
      xB,yB,zB=coordsB[1],coordsB[2],coordsB[3]
      xB,yB,zB=float(xB),float(yB),float(zB)
      qB=float(efragB[j].split()[1])
      r=sqrt((xA-xB)**2+(yA-yB)**2+(zA-zB)**2)
      r*=1e-10 # Angstroms to meters
      # Since E=k*qQ/r and
      # Taking dielectric constant 3 times that of vaccum
      # Coulomb Energy in eV
      myJeV+=qA*qB/(3*4*pi*((8.85418781762e-12)/1.60217646e-19)*r) # in eV
      #print '%s eV' % myJeV
      # Coulomb Energy in Nm (Joules)
      # Since E=k*qQ/r and k=8.99e9 Nm**2C**-2
      qA*=1.602176565e-19 # e to C
      qB*=1.602176565e-19 # e to C
      myJJ+=8.99e9*qA*qB/(3*r) # in Joules
      #print '%s Joules' % myJJ
  return myJeV, myJJ

# now print Coulomb Energy
myJEeV,myJEJ=myCoulomb(cA,cB,eA,eB)      
print '%s eV' % myJEeV
myJEJ*=6.24150974e18 # now in eV
# Now print this only if considerably different
if abs(abs(myJEJ)-abs(myJEeV)) > 0.01:
  print '%s eV (Joule converted)' % myJEJ


