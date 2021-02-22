#! /bin/sh
# Script for modifying the zindo input file generated from the babel. It requires
# xyz file (previously gaussian log files)
# This script is written by Yesudas on August 31st, 2010.
# Modified on Jan 18, 2011. Modified on Feb 11, 2011.
# Modified on Mar 10, 2011 by Laxman Pandey
# k=active space
# nof - Number of Files
# name = name of file
#
script=`echo $0 | sed 's|.*\/\([^\.]*\)|\1|g'`
#
while [ $# -lt 3 ]; do
  echo -e "\nThree Requirements:\n\tA) # of orbitals in CI-space (active space)"
  echo -e "\tB) # of files\n\t\t1: To convert single file\n\t\t2: To convert multiple files"
  echo -e "\tC) choice 1: follows up with 'filename'"
  echo -e "\t   choice 2: follows up with 'all'"
  echo -e "\nUsage: $script 'active space size' '1 or 2' 'filename or all'"
  echo -e "\t(eg. $script 10 1 filename.xyz)\n"
  exit
done
#
k=$1; nof=$2; name=$3
pvn1='1.000000'
#
echo -e "\nParameters given/taken: $1 $2 $3\n"
#
if [ $nof == '2' ] && [ $name == 'all' ]; then
  #for filename in $(ls *.log)
  for filename in $(ls *.opt.xyz); do
    #input=`basename $filename .log`
    input=`basename $filename .opt.xyz`
    #babel -i g03 $filename -o zin $input.inp::q
    babel -i xyz $filename -o zin $input.inp
    nel=`grep 'NEL' $input.inp | awk '{print $4}'`
    nel=`grep 'NEL' $input.inp | awk '{print $4}'`
    homo=`echo "$nel / 2" | bc`
    lumo=`echo "$homo + 1" | bc`
    homo_1=`echo "$homo - 1" | bc`
    homo_2=`echo "$homo - 2" | bc`
    homo_k=`echo "$homo - $k + 1" | bc`
    lumo_1=`echo "$lumo + 1" | bc`
    lumo_k=`echo "$lumo + $k - 1" | bc`
    iprint=`grep  "IPRINT" $input.inp`
    sed -i 's/'"$iprint"'/ IPRINT         -1   ITMAX       150   SCFTOL  0.000001/' $input.inp
    dynal=`grep 'DYNAL' $input.inp`
    p=5000
    dynalmod=`(IFS=' ';set -- $(echo $dynal) ;printf "%9s %1s %5s %4s %4s %4s %4s %4s %4s" $1 $2 $3 $4 $5 $6 $7 $p $lumo_k)`
    sed -i 's/'"$dynal"'/'"$dynalmod"'/' $input.inp
    infa=`grep "INTFA" $input.inp`
    sed -i 's/'"$infa"'/'"`printf "%9s%3s%10s%9s%10s%9s%10s%10s\n" {'INTFA(1)',' = ',$pvn1,'1.267000','0.585000','1.00000',$pvn1,$pvn1}`"'/' $input.inp
    #sed -i 's/'"$infa"'/ INTFA(1) =   1.000000 1.267000  0.585000  1.00000  1.000000  1.000000/' $input.inp
    sed -i 's/zindo/'"$filename"'/' $input.inp:Requirements
    sed -i '$d' $input.inp
    sed -i '$d' $input.inp
    sed -i '$d' $input.inp
    sed -i '$d' $input.inp
    sed -i '$d' $input.inp
    sed -i '$d' $input.inp
    sed -i '$d' $input.inp
    printf "%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s\n" {'5','5','5000','1','0','0','0','1','5000','1','5000'} >> $input.inp
    #echo "    5    5 5000    1    0    0    0    1 5000    1 5000" >> $input.inp
    echo "  -80000.0 0.000000" >> $input.inp
    echo " " >> $input.inp
    printf "%5s %4s %4s\n" 1 $homo $homo >> $input.inp
    printf "%5s %4s %4s %4s %4s\n" 21 $homo_k $homo $lumo $lumo_k >> $input.inp
    printf "%5s %4s %4s %4s %4s\n" 1 $homo_1 $homo_1 $homo $lumo >> $input.inp
    printf "%5s %4s %4s %4s %4s\n" 21 $homo_k $homo $lumo $lumo_k >> $input.inp
    printf "%5s %4s %4s %4s %4s\n" 1 $homo_1 $homo_1 $homo $lumo_1 >> $input.inp
    printf "%5s %4s %4s %4s %4s\n" 21 $homo_k $homo $lumo $lumo_k >> $input.inp
    printf "%5s %4s %4s %4s %4s %4s %4s\n" 1 $homo_2 $homo_2 $homo $homo $homo_1 $lumo >> $input.inp
    printf "%5s %4s %4s %4s %4s\n" 21 $homo_k $homo $lumo $lumo_k >> $input.inp
    printf "%5s %4s %4s %4s %4s\n" 1 $homo_1 $homo_1 $lumo $lumo >> $input.inp
    printf "%5s %4s %4s %4s %4s\n" 21 $homo_k $homo $lumo $lumo_k >> $input.inp
    echo " " >> $input.inp
    echo ' $END' >> $input.inp
  done
elif [ $nof == "1" ]; then
  #echo 'Enter file name' 
  #read filename 
  #input=`basename $name .log`
  filename=$name
  input=`basename $filename .opt.xyz`
  #babel -i g03 $filename -o zin $input.inp
  babel -i xyz $filename -o zin $input.inp
  nel=`grep 'NEL' $input.inp | awk '{print $4}'`
  homo=`echo "$nel / 2" | bc`
  lumo=`echo "$homo + 1" | bc`
  homo_1=`echo "$homo - 1" | bc`
  homo_2=`echo "$homo - 2" | bc`
  homo_k=`echo "$homo - $k + 1" | bc`
  lumo_1=`echo "$lumo + 1" | bc`
  lumo_k=`echo "$lumo + $k - 1" | bc`
  iprint=`grep  "IPRINT" $input.inp`
  sed -i 's/'"$iprint"'/ IPRINT         -1   ITMAX       150   SCFTOL  0.000001/' $input.inp

  #sed -i 's/'"$iprint"'/ IPRINT         -1   ITMAX       150   SCFTOL  0.000001/' $input.inp
  dynal=`grep 'DYNAL' $input.inp`
  p=5000
  dynalmod=`(IFS=' ';set -- $(echo $dynal) ;printf "%9s %1s %5s %4s %4s %4s %4s %4s %4s" $1 $2 $3 $4 $5 $6 $7 $p $lumo_k)`
  sed -i 's/'"$dynal"'/'"$dynalmod"'/' $input.inp
  infa=`grep "INTFA" $input.inp`
  sed -i 's/'"$infa"'/'"`printf "%9s%3s%10s%9s%10s%9s%10s%10s\n" {'INTFA(1)',' = ',$pvn1,'1.267000','0.585000','1.00000',$pvn1,$pvn1}`"'/' $input.inp
  #sed -i 's/'"$infa"'/ INTFA(1) =   1.000000 1.267000  0.585000  1.00000  1.000000  1.000000/' $input.inp
  sed -i 's/zindo/'"$filename"'/' $input.inp
  sed -i '$d' $input.inp 
  sed -i '$d' $input.inp 
  sed -i '$d' $input.inp 
  sed -i '$d' $input.inp 
  sed -i '$d' $input.inp 
  sed -i '$d' $input.inp 
  sed -i '$d' $input.inp 
  #echo "    5    5 5000    1    0    0    0    1 5000    1 5000" >> $input.inp
  printf "%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s%5s\n" {'5','5','5000','1','0','0','0','1','5000','1','5000'} >> $input.inp
  echo "  -80000.0 0.000000" >> $input.inp
  echo " " >> $input.inp
  printf "%5s %4s %4s\n" 1 $homo $homo >> $input.inp
  printf "%5s %4s %4s %4s %4s\n" 21 $homo_k $homo $lumo $lumo_k >> $input.inp
  printf "%5s %4s %4s %4s %4s\n" 1 $homo_1 $homo_1 $homo $lumo >> $input.inp
  printf "%5s %4s %4s %4s %4s\n" 21 $homo_k $homo $lumo $lumo_k >> $input.inp
  printf "%5s %4s %4s %4s %4s\n" 1 $homo_1 $homo_1 $homo $lumo_1 >> $input.inp
  printf "%5s %4s %4s %4s %4s\n" 21 $homo_k $homo $lumo $lumo_k >> $input.inp
  printf "%5s %4s %4s %4s %4s %4s %4s\n" 1 $homo_2 $homo_2 $homo $homo $homo_1 $lumo >> $input.inp 
  printf "%5s %4s %4s %4s %4s\n" 21 $homo_k $homo $lumo $lumo_k >> $input.inp
  printf "%5s %4s %4s %4s %4s\n" 1 $homo_1 $homo_1 $lumo $lumo >> $input.inp
  printf "%5s %4s %4s %4s %4s\n" 21 $homo_k $homo $lumo $lumo_k >> $input.inp
  echo " " >> $input.inp
  echo ' $END' >> $input.inp
fi
echo -e "\nAll Done"
echo " Thank you for using $script"
