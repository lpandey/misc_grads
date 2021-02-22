#!/usr/bin/env python
#Filename: transcomp_gs.py
#Laxman Pandey 6/4/2010

# reads output from gaussian tddft (closed shell) calculation and
# prints frontier MOs contributing to electronic transitions

from __future__ import division
import os, sys
from subprocess import Popen, PIPE 
from itertools import islice

#### ASSIGN APPROPRIATE VALUES ####
au2dby=0.393456 #1dby = 0.393456 a.u.
ev2invcm=8065.73 # 1 eV = 8065.73 cm**-1
fwahmx=0.30 # gaussian full-width-at-half-maximum (eV)
#fwahmx=0.35 # gaussian full-width-at-half-maximum (eV)

def myround(num):
  if (num > 0): return int(num + 0.5)
  else: return int(num - 0.5)

def seek_to_line(f, n):
  for ignored_lines in islice(f, int(n)-1):
    pass #skip beginning n lines

def mymos(anel,exfrm,exto):
  hm,hm_m1,hm_m2,hm_m3,hm_m4,hm_m5,hm_m6,hm_m7,hm_m8,hm_m9,hm_m10, \
  hm_m11,hm_m12,hm_m13,hm_m14,hm_m15,hm_m16,hm_m17,hm_m18,hm_19,hm_m20=\
  anel,anel-1,anel-2,anel-3,anel-4,anel-5,anel-6,anel-7,anel-8,anel-9, \
  anel-10,anel-11,anel-12,anel-13,anel-14,anel-15,anel-16,anel-17,anel-18,\
  anel-19,anel-20

  lm,lm_p1,lm_p2,lm_p3,lm_p4,lm_p5,lm_p6,lm_p7,lm_p8,lm_p9,lm_p10, \
  lm_p11,lm_p12,lm_p13,lm_p14,lm_p15 = \
  anel+1,anel+2,anel+3,anel+4,anel+5,anel+6,anel+7,anel+8,anel+9,\
  anel+10,anel+11,anel+12,anel+13,anel+14,anel+15,anel+16

  if exfrm == hm: exfrm = 'HOMO'
  elif exfrm == hm_m1: exfrm = 'HOMO-1'
  elif exfrm == hm_m2: exfrm = 'HOMO-2'
  elif exfrm == hm_m3: exfrm = 'HOMO-3'
  elif exfrm == hm_m4: exfrm = 'HOMO-4'
  elif exfrm == hm_m5: exfrm = 'HOMO-5'
  elif exfrm == hm_m6: exfrm = 'HOMO-6'
  elif exfrm == hm_m7: exfrm = 'HOMO-7'
  elif exfrm == hm_m8: exfrm = 'HOMO-8'
  elif exfrm == hm_m9: exfrm = 'HOMO-9'
  elif exfrm == hm_m10: exfrm = 'HOMO-10'
  elif exfrm == hm_m11: exfrm = 'HOMO-11'
  elif exfrm == hm_m12: exfrm = 'HOMO-12'
  elif exfrm == hm_m13: exfrm = 'HOMO-13'
  elif exfrm == hm_m14: exfrm = 'HOMO-14'
  elif exfrm == hm_m15: exfrm = 'HOMO-15'
  elif exfrm == hm_m16: exfrm = 'HOMO-16'
  elif exfrm == hm_m17: exfrm = 'HOMO-17'
  elif exfrm == hm_m18: exfrm = 'HOMO-18'
  elif exfrm == hm_m19: exfrm = 'HOMO-19'
  elif exfrm == hm_m20: exfrm = 'HOMO-20'
  if exto == lm: exto = 'LUMO'
  elif exto == lm_p1: exto = 'LUMO+1'
  elif exto == lm_p2: exto = 'LUMO+2'
  elif exto == lm_p3: exto = 'LUMO+3'
  elif exto == lm_p4: exto = 'LUMO+4'
  elif exto == lm_p5: exto = 'LUMO+5'
  elif exto == lm_p6: exto = 'LUMO+6'
  elif exto == lm_p7: exto = 'LUMO+7'
  elif exto == lm_p8: exto = 'LUMO+8'
  elif exto == lm_p9: exto = 'LUMO+9'
  elif exto == lm_p10: exto = 'LUMO+10'
  elif exto == lm_p11: exto = 'LUMO+11'
  elif exto == lm_p12: exto = 'LUMO+12'
  elif exto == lm_p13: exto = 'LUMO+13'
  elif exto == lm_p14: exto = 'LUMO+14'
  elif exto == lm_p15: exto = 'LUMO+15'
  return exfrm, exto

def mymoments(file, tednst, tedlst):
  for line in file:
    line = line.strip()
    if not line: continue
    #print line
    if (line.startswith('Ground') or line.startswith('ground')): break
    teddat = line.split()
    if not (len(teddat) == 5 or len(teddat) == 6): continue
    if not teddat[0].isdigit(): continue
    #print teddat
    assert type(eval(teddat[1])) == float
    assert type(eval(teddat[2])) == float
    assert type(eval(teddat[3])) == float
    assert type(eval(teddat[4])) == float
    if len(teddat) == 6:
      assert type(eval(teddat[5])) == float
      tedos = teddat[5]
    elif len(teddat) == 5:
      tedos = teddat[4]
    #print tedos
    tedst,tedx,tedy,tedz = teddat[0],teddat[1],teddat[2],teddat[3]
    # now transition electric dipole data gathering
    tedlst.append('%s %s %s %s %s' %(tedst, tedx, tedy, tedz, tedos))
  return tedlst

def mytransition(tedstr):
  lst = tedstr.split()
  #print lst
  ftst,rtx,rty,rtz,fosc = lst[0],lst[1],lst[2],lst[3],lst[4]
  ftst,fosc = int(ftst),float(fosc)
  #ftst = int(ftst)
  rtx,rty,rtz = float(rtx),float(rty),float(rtz)
  dip_tot_au = (rtx**2 + rty**2 + rtz**2)**0.5
  ftx,fty,ftz = rtx/au2dby,rty/au2dby,rtz/au2dby
  dip_tot_dby = (ftx**2 + fty**2 + ftz**2)**0.5
  #print ftst, ftx, fty, ftz, dip_tot_dby, fosc      
  return (ftst, ftx, fty, ftz, dip_tot_dby, fosc)      
  #return (ftst, ftx, fty, ftz, dip_tot_dby)      

def my_split(s, seps):
  res = [s]
  for sep in seps:
    s, res = res, []
    for seq in s:
      res += seq.split(sep)
  return res
    
#if len(sys.argv) < 3:
if len(sys.argv) < 2:
  print "\n Grabs important data from gaussian tddft (closed shell) calculation\n"
  print " Usage: transcomp_gs.py filenames"
  print "  (eg. transcomp_gs.py file1.log file2.log)"
  sys.exit()
  #logfs=raw_input(' Pass me tddft filenames:')
  #logfs = logfs.split()
  #plev=raw_input('\n Print level [S1(default), S1-3, S1-5, ST1-3, ST1-5, all] :')
else:
  logfs = sys.argv[1:]; #plev = sys.argv[2]

for logf in logfs:
  myfile=logf.replace('.log','')
  tednst, estnst, crntst = 0, 0, 0
  tedlst, estlst, estlns, dmlst, esoslst = [], [], [], [], []
  allconflst = []
  fort44 = []; mystring = ' %s ' % logf
  #print; #print mystring.center(60, '#') # final print
  #grep strings assuming user of Perl regular expressions -P
  chkstr='"Normal\s+[Tt]ermination"'
  nelstr='"[Aa]lpha\s+[Ee]lectrons"'
  tedstr='"Ground\s+[Tt]o\s+[Ee]xcited\s+[Ss]tate\s+[Tt]ransition\s+[Ee]lectric\s+[Dd]ipole\s+[Mm]oments\s+"'
  eststr='"[Ee]xcited\s+[Ss]tate\s+\d{1,2}:\s+"'
  #eststr='"[Ee]xcited\s+[Ss]tate\s+\d{1,2}:\s+[ST][ir][ni][gp]let"'

  #wt=open('%s.txt' %os.path.basename(sys.argv[0]).replace('.py',''), 'a')
  rlog='%s' % logf
  #cwd=os.getcwd()
  #set up grep commands
  chkcmd='grep -P %s %s' %(chkstr, rlog)
  nelcmd='grep -Pm 1 %s %s' %(nelstr, rlog)
  tedcmd='grep -Pnm 1 %s %s' % (tedstr, rlog)
  estcmd='grep -Pn %s %s' %(eststr, rlog)

  if os.path.exists(rlog):
    #print '%s exists in %s' % (rlog, cwd) 
    #print rlog
    #now check job ended normally
    nrmlend=Popen('%s' %chkcmd, shell=True, stdout=PIPE).stdout
    #nrmlend=nrmlend.read().strip().split()[0]
    nrmlend=nrmlend.read().strip().split()
    #print nrmlend
    if len(nrmlend) < 1:
      print '%s did not terminate normally, or still running' %rlog
      #wt.write('%s DNTN\n' % (logf))
      continue
 
    #now check for closed shell HF or KS
    #nelln=os.system(nelcmd) #different way NA here
    nel=Popen('%s' %nelcmd, shell=True, stdout=PIPE).stdout
    nel=nel.read().split()
    #print nel
    anel, bnel = eval(nel[0]), eval(nel[3])
    #wt.write('%s %s %s\n' %(logf, anel, bnel))
    if anel != bnel:
      print rlog, '# alpha electrons != # beta electrons'
      #wt.write('%s %s != %s\n' %(logf, anel, bnel))
      #sys.exit()
      continue
      
    #now Transition electric dipole moments
    tedln=Popen('%s' %tedcmd, shell=True, stdout=PIPE).stdout
    tedln=tedln.read().split()[0][:-1] #get rid of :
    tedln=int(tedln)
    #print tedln
    rf=open('%s' %rlog, 'r') #open log file rlog
    #seek_to_line(rf, tedln)
    seek_to_line(rf, tedln + 1)
    mymoments(rf, tednst, tedlst)
    rf.close()
    #print tedlst # print moments list
    #print "\n%-7s%7s%12s%12s%12s%10s" %('State','tdpx','tdpy','tdpz','tdp(tot)','Osc') # final print
    
    
    for i in range(len(tedlst)):
      tedstr = tedlst[i]
      #print tedstr
      #mytransition(tedstr)
      tst,tx,ty,tz,tt,osc = mytransition(tedstr)
      #tst,tx,ty,tz,tt = mytransition(tedstr)
      #print "%-2d%12.4f%12.4f%12.4f%12.4f%10.4f" %(tst,tx,ty,tz,tt,osc) # final print
      dmlst.append("%-2d%12.4f%12.4f%12.4f%12.4f%10.4f" %(tst,tx,ty,tz,tt,osc))
 
    #for i in range(1,len(tedlst)):
    #  data['state%s' %i] = {}
      
    #print #final print
    #read excited state energies and configurations
    estln=Popen('%s' %estcmd, shell=True, stdout=PIPE).stdout
    estln=estln.readlines()
    for i in range(len(estln)):
      estlns.append(estln[i].split()[0][:-1])
    #print len(estlns), estlns
    #print 'estlns', estlns
 
    for i in range(len(estlns))[:1]:
      #print 'I am here?'
      rf=open('%s' %rlog, 'r') #open log file rlog
      #print int(estlns[i])
      seek_to_line(rf, estlns[i])
      for esline in rf:
        #print esline
        esline=esline.strip().replace('->','');
        esline=esline.strip().replace('<-','');
        esline=esline.strip().replace('A','');
        esline=esline.strip().replace('B','');
        if not esline: continue
        #print esline
        esendat=esline.split();
        if esendat[0].startswith('******'): break
        #print 'len esendat', len(esendat)
        if not (len(esendat)==9 or len(esendat)==10 or len(esendat)==3):
          continue
        if not (esendat[2][:-1].isdigit() or esendat[1].isdigit()):
          continue
        #print esendat
        try:
          esst=esendat[2][:-1]; estype=esendat[3]; esev=esendat[4];
          esnm=esendat[6]; esos=esendat[8][2:]
          esst = int(esst); esev = float(esev)
          esnm = float(esnm); esos = float(esos); esinvcm=esev*ev2invcm/1000.0
          #print 'esos', esos
          #print "%-3d%7s%7.3f%5d%7.2f" %(esst,estype[:-2],esev,esnm,esos) 
          #print "%-3d%12s%7.3f%6.1f%5d%7.2f" %(esst,estype,esev,esinvcm,esnm,esos) #final print 
          esoslst.append("%-3d%12s%7.3f%6.1f%5d%7.2f" %(esst,estype,esev,esinvcm,esnm,esos))
          #esoslst.append("%-7s%3d%7.3f%7d%6.1f%7.2f" %(estype[:-2],esst,esev,esinvcm,esnm,esos))
          #print "I am here"
          fort44.append("%7.3f%6.2f%6.1f" %(esev, esos, esinvcm))
          #print fort44 
          crntst=esst
        except: pass
        #except: continue 
        #print 'esendat', esendat

        if (len(esendat) == 3 and esendat[0].isdigit() and esendat[1].isdigit()):
        #if (len(esendat) == 3 and esendat[0][:-1].isdigit() and esendat[1][:-1].isdigit()):
          #print 'esendat', esendat
          try:
            type(eval(esendat[2])) == float
            #print 'esendat', esendat
            exfrm=esendat[0];exto=esendat[1];excof=esendat[2]
            exfrm,exto,pexcof=int(exfrm),int(exto),2*float(excof)**2*100
            #print '  %-8s%-8s%3d' %(exfrm, exto, pexcof)
            exfrm, exto = mymos(anel,exfrm,exto)
            #print '  %-8s%-8s%9s' %(exfrm, exto, excof)
            #print '   %-8s%-8s%3d' %(exfrm, exto, pexcof)
            #print '   %-s-->%s(%s)' %(exfrm, exto, myround(pexcof)) #final print
            #esoslst.append('%-8s%-8s%9s'%(exfrm,exto,excof))
            #mylist.append('%-8s%-8s%3d' %(exfrm, exto, excof))
            allconflst.append("%s %-s-->%s(%s)" %(crntst, exfrm, exto, myround(pexcof)))
          except: continue
      rf.close()
    #print fort44
    if len(fort44):
      fort45=open('%s.dat44' % rlog.replace('.log',''), 'w')
      fort45.write("%3d\n" % len(fort44))
      minE=eval(fort44[0].split()[0])-1
      if minE < 0: minE = 0.05
      maxE=eval(fort44[len(fort44)-1].split()[0])+1
      fort45.write("%6.2f%7.3f%7.3f\n" % (fwahmx, minE, maxE))
      for i in range(len(fort44)):
        fort45.write("%7.3f%6.2f%7d%7d\n" %( eval(fort44[i].split()[0]),
                      eval(fort44[i].split()[1]), 1240.0/eval(fort44[i].split()[0]),
                      eval(fort44[i].split()[2])*1000))    
        #fort45.write("%s\n" % fort44[i])
      fort45.close()
    #print dmlst 
    #print esoslst
    #print allconflst
    stateconflsts = []
    for i in range(len(dmlst)):
      stateconflst = []
      for j in range(len(allconflst)):
        if int(allconflst[j].split()[0]) == int(dmlst[i].split()[0]):
          stateconflst.append(allconflst[j].split()[1])
      stateconflsts.append(stateconflst)
    #print 'dmlst', dmlst
    #print 'esoslst', esoslst
    #print 'stateconflsts', stateconflsts  
    assert len(dmlst) == len(esoslst) == len(stateconflsts)
    #wf=open('%s.dat' %myfile, 'w')
    #wf.write("%-20s" %myfile)
    print myfile,
    for k in range(len(dmlst)):
    #for k in range(1):
      state =int(dmlst[k].split()[0])
      mu    =float(dmlst[k].split()[4])
      f     =float(dmlst[k].split()[5])
      en    =float(esoslst[k].split()[2]) 
      mult  =esoslst[k].split()[1]      
      #print " %2s %s %6.3f %5.2f %5.2f" %(state,mult,en,mu,f,conf)
      #print "%-6.3f %5.2f %5.2f" %(en,mu,f)
      for ci in range(len(stateconflsts[k])):
        conf=stateconflsts[k][ci]
        #if "HOMO-->LUMO(" in conf:
          #print " %2s %s %5.2f %5.2f %5.2f %s" %(state,mult,en,mu,f,conf)
        print " %2s %s %6.3f %5.2f %5.2f %s" %(state,mult,en,mu,f,conf)
        #wf.write(" %2s %s %5.2f %5.2f %5.2f %s\n" %(state,mult,en,mu,f,conf))
    #wf.close()

  else:
    print "\n %s not in current directory !!!\n" %rlog

print "\nGaussian FW@HM used = %g eV\n" % fwahmx
