#!/bin/sh
# Laxman Pandey
# Check gaussian optimization status
while [ $# -lt 1 ]; do
    #scriptname=`echo $0 | sed 's|.*\/\([^\.]*\)\(\..*\)$|\1|g'`
    scriptname=`echo $0 | sed 's|.*\/\([^\.]*\)|\1|g'`
    echo -e "\nUsage: $scriptname logfile.log\n"
    exit
done

file=$1 
egrep "out of|SCF Don|Converged| NO | YES|exceeded| Optimization " $file | grep -v '\\\\'
