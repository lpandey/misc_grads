#!/bin/env python
#written for pbs environment on FORCE cluster
#g03 or g09 versions
#gjf2pbs.py by Laxman Pandey

import os, sys, glob

if len(sys.argv) < 8:
  print "\n ******** Creates submit scripts for gaussian input files ********"
  print " Usage: gjf2pbs.py queue gsver #nodes #procs #pmem walltimehrs filenames"
  print " (eg. gjf2pbs.py force g09 1 2 2000 100 test1.gjf test2.gjf\n"
  sys.exit()

gsver,nnodes,nprocs,pmem=sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5]
queue,walltimehrs = sys.argv[1],sys.argv[6]
assert (gsver=='g03' or gsver=='g09'), 'gaussian version error: (Use: g03 or g09)'
assert (queue=='force' or queue=='optimus'), 'queue error: (Use: force or optimus)'
try:
  eval(nnodes); eval(nprocs); eval(pmem); eval(walltimehrs)
except:
  print "Error in #nodes #procs #pmem walltimehrs (eg. 1 2 1500 500)"
  sys.exit()

#for file in glob.glob('*.gjf'):   
for file in sys.argv[7:]:
  gjff=file.replace('.gjf','')
  wf2=open('%s.pbs' %gjff, 'w')
  wf2.write('#!/bin/bash\n')
  wf2.write('# This file is prepared by gjf2pbs.py\n\n')
  wf2.write('## Specify which queue to send to\n')
  wf2.write('#PBS -q %s\n\n' % queue)
  wf2.write('## Merge standard error and standard output\n')
  wf2.write('#PBS -j oe\n\n')
  wf2.write('## %s node(s) with %s proc(s) per node\n' %(nnodes,nprocs))
  wf2.write('#PBS -l nodes=%s:ppn=%s\n\n' %(nnodes,nprocs))
  wf2.write('## specifying %s mb memory per procs\n' %pmem)
  wf2.write('#PBS -l pmem=%smb\n\n' %pmem)
  wf2.write('## calcualtion name\n#PBS -N %s\n\n' %gjff)
  wf2.write('## specify %s hours CPU walltime\n' %walltimehrs)
  wf2.write('#PBS -l walltime=%s:00:00\n\n' %walltimehrs)
  wf2.write('## Call PBS to send email when job ends\n')
  wf2.write('#PBS -M lpandey@gmail.com\n#PBS -m ae\n\n')
  wf2.write('## export environment defined at submission to compute node\n')
  wf2.write('#PBS -V\n\n## source the Gaussian environment\n')
  if gsver == 'g03':
    wf2.write('export %sroot=/usr/local/packages/gaussian/G03E01\n' %gsver)
  elif gsver == 'g09':
    wf2.write('export %sroot=/usr/local/packages/gaussian/G09B01LINDA\n' %gsver)
  wf2.write('. ${%sroot}/%s/bsd/%s.profile\n' %(gsver,gsver,gsver))
  wf2.write('export PATH=$TMPDIR:$PATH\n\n')
  wf2.write('## Change to working directory\ncd $PBS_O_WORKDIR\n\n')
  wf2.write('## scratch dir name using Job name and Job ID from PBS\n')
  #wf2.write('if [ `mount | grep /tmp` ]; then\n  if [ ! -d /tmp/$USER ]; then\n')
  #wf2.write('    mkdir /tmp/$USER\n  fi\n')
  #wf2.write('  myscratch=/tmp/$USER/${PBS_JOBNAME}_$(echo $PBS_JOBID | cut -d"." -f1);\n')
  #wf2.write('else\n  ')
  #wf2.write('myscratch=/nv/hp18/$USER/scratch/${PBS_JOBNAME}_$(echo $PBS_JOBID | cut -d"." -f1);\n')
  #wf2.write('fi\n')
  wf2.write('myscratch=/nv/hp18/$USER/scratch/${PBS_JOBNAME}_$(echo $PBS_JOBID | cut -d"." -f1);\n')
  wf2.write('export myscratch\n\n')
  wf2.write('## Creating scratch directoy\n')
  wf2.write('mkdir $myscratch\nexport GAUSS_SCRDIR=$myscratch\n\n')
  wf2.write('/bin/rm -f $PBS_JOBNAME.node\n')
  wf2.write('for host in $(cat $PBS_NODEFILE); do\n')
  wf2.write('  echo "${host}" >> $PBS_JOBNAME.node\ndone\n\n')
  wf2.write('## set up the Gaussian running environment\n')
  wf2.write('export NODES=\\"`cat $PBS_JOBNAME.node`\\"\n')
  wf2.write('export GAUSS_LFLAGS="-nodefile $PBS_JOBNAME.node -vv"\n\n')
  wf2.write('## Run the job in the scratch directories\n')
  wf2.write('cat $PBS_JOBNAME.node\nexport GAUSS_SCRDIR=$myscratch\n')
  wf2.write('date\n%s $PBS_JOBNAME.gjf\ndate\n\n' %gsver)
  wf2.write('## Clean up after by removing the scratch directories\n')
  wf2.write('/bin/rm -r -f $myscratch\n')
  print ' %s.pbs written' %gjff

