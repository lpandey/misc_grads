#!/usr/bin/env python
# Laxman Pandey 07/06/2011

# set-up omegaTuning jobs (Long Range Corrected Density Functionals)

import sys, os
from subprocess import Popen, PIPE
from itertools import islice

if len(sys.argv) != 3:
  print "\n  Usage: omegaTuning_g09_DA_data.py filename.log NumberOfJobs\n"
  sys.exit(1)

logf = sys.argv[1]
numjobs = int(sys.argv[2]) 

# Define grep execution function
def myGrep(cmd):
  grepOut = Popen(cmd, shell=True, stdout=PIPE).stdout
  return grepOut

# Define grep command function
def myGrepCmd(cmdStr, file, flag=''):
  grepCmd = 'grep %s "%s" %s' % (flag, cmdStr, file)
  return grepCmd 

# Define list print function
def printList(myList, newline='yes'):
  for i in range(len(myList)):
    if newline == 'yes': print myList[i]
    elif newline == 'no': print myList[i], 

# Define skip to a specific line in a file function
def seek_to_line(f, n):
  for ignored_lines in islice(f, int(n)-1):
    pass # skip beginning n lines

# Define function to read and extract eigenvalues from lines
def myOrbitals(file, aOcclst, aVirtlst, bOcclst, bVirtlst):
  for line in file:
    line = line.strip('\n').strip()
    if not line: continue  # skip out empty lines
    if 'Condensed to atoms' in line: break # break out of loop
    orbData = line.split()
    try: type(eval(orbData[4])) == float # check eigenvalues lines
    except: continue # skip otherwise
    if (orbData[0] == 'Alpha' and orbData[1] == 'occ.'):
      alphaOcc = orbData[4:] # only take the eigenvalues
      aOcclst.extend(alphaOcc) # append each item to list
    elif (orbData[0] == 'Alpha' and orbData[1] == 'virt.'):
      alphaVirt = orbData[4:]
      aVirtlst.extend(alphaVirt)
    elif orbData[0] == 'Beta' and orbData[1] == 'occ.':
      betaOcc = orbData[4:]
      bOcclst.extend(betaOcc)
    elif orbData[0] == 'Beta' and orbData[1] == 'virt.':
      betaVirt = orbData[4:]
      bVirtlst.extend(betaVirt)
  return aOcclst, aVritlst, bOcclst, bVirtlst

# Define HOMO and LUMO determination function
def getEigenvalues(aOcc, aVirt, bOcc, bVirt):
  try:
    if len(bOcc) > 0 and len(bVirt) > 0:
      aHOMO, aLUMO, bHOMO, bLUMO = float(aOcc[-1]), float(aVirt[0]), \
                                   float(bOcc[-1]), float(bVirt[0])
      if aHOMO >= bHOMO: HOMO = aHOMO
      else: HOMO = bHOMO
      if aLUMO  <= bLUMO: LUMO = aLUMO
      else: LUMO = bLUMO
    else:
      HOMO, LUMO = float(aOcc[-1]), float(aVirt[0])
  except:
    pass
  return HOMO, LUMO

# grep command strings
jobsTotStr = '[Ii]nitial\s+[Cc]ommand:'
normEndStr = '[Nn]ormal\s+[Tt]ermination'
scfStr     = 'SCF Done'
popLnStr='[Pp]opulation\s+[Aa]nalysis\s+[Uu]sing\s+[Tt]he\s+SCF\s+[Dd]ensity'

# grep commands
jobsTotCmd = myGrepCmd(jobsTotStr, logf, flag='-Pc') # count
normEndCmd = myGrepCmd(normEndStr, logf, flag='-Pn') # line num
scfEnCmd   = myGrepCmd(scfStr, logf)
popLnCmd   = myGrepCmd(popLnStr, logf, flag='-Pn')   # line num

normEndLnlst, scflst, popLnlst = [], [], []
HOMOlst, LUMOlst = [], []

# Check normal termination of jobs
jobsTot = myGrep(jobsTotCmd).read().strip('\n') # use myGrep
#print "\n%s" % jobsTot, "linked jobs found..."

normEndLns = myGrep(normEndCmd).readlines() # use myGrep 

for i in range(len(normEndLns)):
  normEndLnlst.append(normEndLns[i].split()[0][:-1]) # cut :

#print " ...with %s Normal ends at lines:" % len(normEndLnlst),
#printList(normEndLnlst, newline='no') # use printList

if len(normEndLnlst) != numjobs:
  print "  Not all jobs ended Normally %s" % len(normEndLnlst)
  exit(1)

#print
# Get SCF energies and SCF calc types
scfEns = myGrep(scfEnCmd).readlines() # use myGrep
for scfEn in range(len(scfEns)):
  scf = scfEns[scfEn].replace('=', '') # remove equals sign
  #scf = scf.split()[2:4] # keeps calc type info
  scf = scf.split()[3]
  #print scf
  scflst.append(scf)

#printList(scflst, newline='yes') # use printList

# Get population information
popLns = myGrep(popLnCmd).readlines() # use myGrep
for i in range(len(popLns)):
  popLnlst.append(popLns[i].split()[0][:-1]) # cut :

#print "\n\n%s sets of Orbital eigenvalues found...\n ...starting at lines:" \
#          % len(popLnlst),
#printList(popLnlst, newline='no') # use printList

# Get population values
for i in range(len(popLnlst)):
  aOcclst, aVritlst, bOcclst, bVirtlst = [], [], [], []
  readf = open(logf, 'r') # open log file to read
  seek_to_line(readf, int(popLnlst[i])) # use seek_to_line
  aOcc, aVirt, bOcc, bVirt  = myOrbitals(readf, aOcclst, aVritlst, 
                              bOcclst, bVirtlst) # use myOrbitals
  readf.close()
  HOMO, LUMO = getEigenvalues(aOcc, aVirt, bOcc, bVirt) # use getEigenvalues
  HOMOlst.append(HOMO)
  LUMOlst.append(LUMO)
  #print '%10.5f %10.5f' % (HOMO, LUMO)  

# Now use collected data for formatted print
mydata = {}
if numjobs == 6:
  jobs = ['Nut:bas1', 'Nut:bas2', 'Cat:bas1', 'Cat:bas2', 'An:bas1', 'An:bas2']
elif numjobs == 3:
  jobs = ['Nut:bas1', 'Cat:bas1', 'An:bas1']

assert len(scflst) == len(HOMOlst) == len(LUMOlst)
#print "\n\n%10s%10s%10s%10s%10s%10s" %('Nut:bas1', 'Nut:bas2', 'Cat:bas1', 'Cat:bas2', 
#      'An:bas1', 'An:bas2\n')
for i in range(len(scflst)):
  mydata[jobs[i]] = (scflst[i], HOMOlst[i], LUMOlst[i])

#print scflst, '\n'
#print HOMOlst, '\n'
#print LUMOlst, '\n'
#print "\n\n", mydata

#sys.exit()

if numjobs == 6:
  print "%18.8f%10.5f%10.5f%18.8f%10.5f%10.5f%18.8f%10.5f%10.5f%18.8f%10.5f%10.5f%18.8f%10.5f%10.5f%18.8f%10.5f%10.5f" % (float(mydata['Nut:bas1'][0]), float(mydata['Nut:bas1'][1]),
        float(mydata['Nut:bas1'][2]), float(mydata['Cat:bas1'][0]),
        float(mydata['Cat:bas1'][1]), float(mydata['Cat:bas1'][2]),
        float(mydata['An:bas1'][0]), float(mydata['An:bas1'][1]),
        float(mydata['An:bas1'][2]), float(mydata['Nut:bas2'][0]),
        float(mydata['Nut:bas2'][1]), float(mydata['Nut:bas2'][2]),
        float(mydata['Cat:bas2'][0]),  float(mydata['Cat:bas2'][1]),
        float(mydata['Cat:bas2'][2]),  float(mydata['An:bas2'][0]),
        float(mydata['An:bas2'][1]),  float(mydata['An:bas2'][2]))
elif numjobs == 3:
  print "%18.8f%10.5f%10.5f%18.8f%10.5f%10.5f%18.8f%10.5f%10.5f" % (
        float(mydata['Nut:bas1'][0]), float(mydata['Nut:bas1'][1]),
        float(mydata['Nut:bas1'][2]), float(mydata['Cat:bas1'][0]),
        float(mydata['Cat:bas1'][1]), float(mydata['Cat:bas1'][2]),
        float(mydata['An:bas1'][0]),  float(mydata['An:bas1'][1]),
        float(mydata['An:bas1'][2]))
#elif numjobs == 1:
#  print "%18.8f%10.5f%10.5f%10.5f%10.5f" % (
#        float(mydata['Nut:bas1'][0]), float(mydata['Nut:bas1'][1]),
#        float(mydata['Nut:bas1'][2]), float(mydata['Nut:bas1'][3]))

#print mydata['Nut:bas1'] 


