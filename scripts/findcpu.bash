#!/bin/bash
# @(#) lists cpu usage of a selection of computers
#
if [ $# -le 0 ]; then
   echo Usage: $0 machinenames
else
for x in $@; do 
   echo 'checking cpu usage on machine' $x;
   ssh $x $PWD/checkcpu.pl;
done
fi
