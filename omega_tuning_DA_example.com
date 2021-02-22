#######################################
Template input (moleculeName_myOmega) file:
#######################################
%nprocs=1
%mem=1000mb
%chk=moleculeName_myOmega.chk
# wB97/6-31G** nosymm Int=ultrafine
# IOp(3/107=0myOmega000000)IOp(3/108=0myOmega000000)
# scf(tight,maxcycles=220)

omega tuning neutral 6-31G** (0.myOmega bohr**-1)

0 1
Coordinates of the molecules here

--Link1--
%nprocs=1
%mem=1000mb
%chk=moleculeName_myOmega.chk
# wB97/6-31G** nosymm Int=ultrafine
# IOp(3/107=0myOmega000000)IOp(3/108=0myOmega000000)
# scf(tight,maxcycles=220) guess=read geom=check

omega tuning cation 6-31G** (0.myOmega bohr**-1)

1 2

--Link1--
%nprocs=1
%mem=1000mb
%chk=moleculeName_myOmega.chk
# wB97/6-31G** nosymm Int=ultrafine
# IOp(3/107=0myOmega000000)IOp(3/108=0myOmega000000)
# scf(tight,maxcycles=220) geom=check

omega tuning anion 6-31G** (0.myOmega bohr**-1)

-1 2

################################################
script to prepare gjf (once one has the template input file).
################################################

for i in $(seq -w 125 5 160);
do
  myname=moleculeName
  cp ${myname}_myOmega.gjf ${myname}_$i.gjf;
  sed -i "s/myOmega/$i/g" ${myname}_$i.gjf;
done

#####################################################
script to extract the values from gaussian output (attached herein)
#####################################################

