#!/bin/sh
#
# @(#) renames output files - D. Price 26/4/04
#
if [ "$#" -ne 2 ]; then
   echo Usage: $0 runname newname
else
   runname=$1;
   newname=$2;
   if [ -f $runname.ev ]; then echo $newname.ev; mv $runname.ev $newname.ev;
   else echo no file $runname.ev; fi
   if [ -f $runname.dat ]; then echo mv $runname.dat $newname.dat; 
   else echo no file $runname.dat; fi
   if [ -f $runname.log ]; then echo mv $runname.log $newname.log;
   else echo no file $runname.log; fi
   if [ -f $runname.sge ]; then echo mv $runname.sge $newname.sge;
   else echo no file $runname.sge; fi
   if [ -f $runname.output ]; then echo mv $runname.output $newname.output;
   else echo no file $runname.output; fi
   if [ -f $runname.errors ]; then echo mv $runname.errors $newname.errors;
   else echo no file $runname.errors; fi
fi
