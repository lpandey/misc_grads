#! /usr/bin/env python
# Laxman Pandey 1.17.2011

# reads Gaussian geometry optimization log file 
# and prints out the HOMO-4 ... LUMO+4 eigenvalues 
import os, sys
from subprocess import Popen, PIPE
from itertools import islice

if len(sys.argv) < 3:
  print "\n Usage: grab_gsopt_HL.py logfile printLevel"
  print " (eg. grab_gsopt_HL.py file.log H-1L+1)\n" 
  #logf = raw_input('\nPass me logfile : ')
  #printlev = raw_input('\nPrint level [HL (default), H-1L+1, H-2L+2, ... H-4L+4] : ')   
  sys.exit(1)

logf = sys.argv[1]; printlev = sys.argv[2]

#print
necmd = 'grep -m 1 %s %s' %('"[Aa]lpha [Ee]lectrons"',logf)
chkcmd = 'grep -c "%s" %s' %('[Oo]ptimization completed', logf)
strtlncmd = 'grep -n "%s" %s' %('Population analysis using the SCF density', logf)
atmcmd = 'grep -m 1 NAtoms= %s' % logf
aoccens, avirens = [], []
boccens, bvirens = [], []

#print "HOMO-4 HOMO-3 HOMO-2 HOMO-1 HOMO LUMO LUMO+1 LUMO+2 LUMO+3 LUMO+4"

def seek_to_line(f, n):
  for ignored_lines in islice(f, n):
    pass # skip beginning n lines
  
if os.path.exists(logf):
  #print logf
  optconfirm=Popen('%s' %chkcmd, shell=True, stdout=PIPE).stdout
  optconfirm=optconfirm.read().strip()
  if int(optconfirm) == 0:
      print '%s Oops, geomtery did not optimize !!!\n'% logf; sys.exit()
  ne = Popen('%s' % necmd, shell=True, stdout=PIPE).stdout
  ne = ne.read().strip().split()
  nalphae, nbetae = int(ne[0]), int(ne[3])
  #print nalphae, nbetae
  if nalphae != nbetae:
      #print "%s !! Unequal Alpha & Beta electrons !!" % logf
      #sys.exit()
      pass
  strtln = Popen('%s' % strtlncmd, shell=True, stdout=PIPE).stdout.read()
  strtln = strtln.strip().split('\n')
  #print strtln
  strtln = strtln[1].split()[0][:-1] # get to second population paragraph
  strtln = int(strtln)
  #print strtln
  rf = open(logf, 'r')
  seek_to_line(rf, strtln)
  # Now find and print the eigenvalues        
  for line in rf:
    line = line.strip('\n').strip()
    if not line: continue
    check = line.split()
    if len(check) < 5 or len(check) > 9: continue # to skip other lines
    try: i = float(check[4])
    except ValueError:
      continue    # skip non-essential lines
    try:
      if check[0] == "Alpha" and len(avirens) <= 5:
        #print 'beta check', check
        if check[1] == 'occ.':
          aoccl = check[4:]
          for n in range(len(aoccl)):
            aoccens.append(aoccl[n])
        if check[1] == 'virt.':
          avirl = check[4:]
          for n in range(len(avirl)):
            avirens.append(avirl[n])
      if check[0] == "Beta" and len(bvirens) <=5:
        #print 'beta check', check
        if check[1] == 'occ.':
          boccl = check[4:]
          for n in range(len(boccl)):
            boccens.append(boccl[n])
        if check[1] == 'virt.':
          bvirl = check[4:]
          for n in range(len(bvirl)):
            bvirens.append(bvirl[n])
    except: continue
    if bvirens:
      if len(bvirens) > 5:
        break
  rf.close()
else:
  print '%s DOES NOT EXIST !!!' %logf
  sys.exit(1)

#print 'aoccens', aoccens
#print 'avirens', avirens
#print 'boccens', boccens
#print 'bvirens', bvirens


# Now assign eigenvalue numbers
aHOMO_m4,aHOMO_m3,aHOMO_m2,aHOMO_m1,aHOMO = \
    float(aoccens[nalphae-5]), float(aoccens[nalphae-4]), \
    float(aoccens[nalphae-3]), float(aoccens[nalphae-2]), \
    float(aoccens[nalphae-1])  

aLUMO,aLUMO_p1,aLUMO_p2,aLUMO_p3,aLUMO_p4 = \
    float(avirens[0]), float(avirens[1]), float(avirens[2]), \
    float(avirens[3]), float(avirens[4])
 
if boccens and bvirens:
  bHOMO_m4,bHOMO_m3,bHOMO_m2,bHOMO_m1,bHOMO = \
    float(boccens[nbetae-5]), float(boccens[nbetae-4]), \
    float(boccens[nbetae-3]), float(boccens[nbetae-2]), \
    float(boccens[nbetae-1])

  bLUMO,bLUMO_p1,bLUMO_p2,bLUMO_p3,bLUMO_p4 = float(bvirens[0]), \
      float(bvirens[1]), float(bvirens[2]), float(bvirens[3]), \
      float(bvirens[4])

# Now evaluate eigenvalues (from alpha and beta)
if boccens and bvirens:
  HOMOs = sorted([aHOMO_m4, aHOMO_m3, aHOMO_m2, aHOMO_m1, aHOMO, \
      bHOMO_m4, bHOMO_m3, bHOMO_m2, bHOMO_m1, bHOMO])[5:]  
  LUMOs = sorted([aLUMO_p4, aLUMO_p3, aLUMO_p2, aLUMO_p1, aLUMO, \
      bLUMO_p4, bLUMO_p3, bLUMO_p2, bLUMO_p1, bLUMO])[:5]  

else:
  HOMOs = [aHOMO_m4, aHOMO_m3, aHOMO_m2, aHOMO_m1, aHOMO] 
  LUMOs = [aLUMO, aLUMO_p1, aLUMO_p2, aLUMO_p3, aLUMO_p4] 

# Now print eigenvalues
#print "H-4 H-3 H-2 H-1 H L L+1 L+2 L+3 L+4"
if printlev == 'H-4L+4':
  print "%s%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f" \
        % (logf, HOMOs[0], HOMOs[1], HOMOs[2], HOMOs[3], HOMOs[4], \
            LUMOs[0], LUMOs[1], LUMOs[2], LUMOs[3], LUMOs[4])
elif printlev == 'H-3L+3':
  print "%s%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f" \
        % (logf, HOMOs[1], HOMOs[2], HOMOs[3], HOMOs[4], \
            LUMOs[0], LUMOs[1], LUMOs[2], LUMOs[3])
elif printlev == 'H-2L+2':
  print "%s%11.5f%11.5f%11.5f%11.5f%11.5f%11.5f" \
        % (logf, HOMOs[2], HOMOs[3], HOMOs[4], LUMOs[0], \
            LUMOs[1], LUMOs[2])
elif printlev == 'H-1L+1':
  print "%s%11.5f%11.5f%11.5f%11.5f" \
        % (logf, HOMOs[3], HOMOs[4], LUMOs[0], LUMOs[1])
else:
  print "%s %s %s" % (logf, HOMOs[4], LUMOs[0])


