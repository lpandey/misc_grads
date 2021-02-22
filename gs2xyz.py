#! /usr/bin/env python2.7
# Laxman Pandey 9.10.2010

# grabs xyz coordinates if geometry optimization (gaussian) ran successfully

import os, sys, re
from subprocess import Popen, PIPE
from itertools import islice

def seek_to_line(f, n):
  for ignored_lines in islice(f, n+1):
    pass # skip beginning n+1 lines

if len(sys.argv) < 2:
  print "\nUsage: gs2xyz.py optimizationlogfile"
  print "  (eg. gs2xyz.py file1.log file2.log)\n"
  sys.exit()

for logf in sys.argv[1:]:
  xyzf = '%s' % logf.replace('.log', '')
  chkcmd = 'grep -c "%s" %s' %('[Oo]ptimization completed', logf)
  strtlncmd = 'grep -n "%s" %s' %('[Oo]ptimization completed', logf)
  atmcmd = 'grep -m 1 NAtoms= %s' % logf
  atoms = 0

  if os.path.exists(logf):
    print logf
    optconfirm=Popen('%s' %chkcmd, shell=True, stdout=PIPE).stdout
    optconfirm=optconfirm.read().strip()
    if int(optconfirm) == 0:
      print ' Oops, geomtery did not optimize\n'
      #sys.exit()
      continue
    natomrf = Popen('%s' % atmcmd, shell=True, stdout=PIPE).stdout
    natomrf = natomrf.read().strip().split()[1]
    #print ' %s atoms found' % natomrf
    natomrf = int(natomrf)
    strtln = Popen('%s' % strtlncmd, shell=True, stdout=PIPE).stdout.read()
    strtln = strtln.strip().split('\n')
    strtln = strtln[0].split()[0][:-1]
    strtln = int(strtln)
    #Now extracting the coordinates
    rf = open(logf, 'r')
    seek_to_line(rf, strtln)
    # Now find and print the coord lines
    wf = open('%s.opt.xyz' %xyzf, 'w')
    wf.write('%d\n\n' % natomrf)
    for line in rf:
      #lines=re.findall(r's+\d{1,3}s+\d{1,3}\s+\d{1}\s+-?\d{1,3}\.\d{1,6}\s+-?\d{1,3}\.\d{1,6}\s+-?\d{1,3}\.\d{1,6}', line)
      line = line.strip('\n').strip()
      if not line: continue
      check = line.split()
      if len(check) != 6: continue
      check = check[0]
      #print check, type(check)
      #if line.startswith('!'):continue
      if not check.isdigit(): continue
      #print line
      atoms = atoms + 1
      line = line.split()
      atom_num,elem,x,y,z=line[0],line[1],line[3],line[4],line[5]
      #print type(elem)
      if elem == '1': elem = 'H'
      if elem == '6': elem = 'C'
      if elem == '7': elem = 'N'
      if elem == '8': elem = 'O'
      if elem == '9': elem = 'F'
      if elem == '14': elem = 'Si'
      if elem == '16': elem = 'S'
      if elem == '17': elem = 'Cl'
      if elem == '34': elem = 'Se'
      #print '%-3s %15s %12s %12s' % (elem, x, y, z)
      wf.write(' %-3s %15s %12s %12s\n' % (elem, x, y, z))
      if atoms == natomrf :
        break
      else:
        pass
    if atoms != natomrf or int(atom_num) != natomrf:
      print ' !!! NUMBER OF ATOMS MISMATCH !!!\n'
      wf.close()
      print '     Deleting opened write file %s.opt.xyz\n' %xyzf
      os.system('rm %s.opt.xyz' %xyzf)
      sys.exit()
    print ' %d atoms written to %s.opt.xyz' % (natomrf,xyzf)
    wf.write('\n')
    rf.close()
    wf.close()
  else:
    print '%s does not exist' %logf
  print

