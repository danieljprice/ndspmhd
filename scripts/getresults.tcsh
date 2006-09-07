#!/bin/tcsh
# Collates results of getav and hsoftvalues into results file
#
if $# != 1 then
   echo Usage: $0 directory
else
   set run=$1
   cd $run
   set date=`date +'_%a_%b_%d_%H:%M:%S'`
   if ( -e results ) then
      echo 'renaming old results -> results'$date
      mv results results$date
   endif
   echo 'using getav...'
   if ( -e crap ) rm crap
   ~/ndspmhd/plot/scripts/getav.pl $run?.output $run??.output > crap
   cat crap
   cat crap | cut -d':' -f 3 > results
   rm crap
   if ( -e hsoft_$run ) then
      cut -d' ' -f 2 hsoft_$run > hsoftvalues
      if ( -e results_$run.txt ) mv results_$run.txt results_$run$date.txt
      paste hsoftvalues results > results_$run.txt
      cp results_$run.txt ../
   else 
      if ( -e hfact_$run ) then
         cut -d' ' -f 3 hfact_$run > hneighvalues
         if ( -e results_$run.txt ) mv results_$run.txt results_$run$date.txt
         paste hneighvalues results > results_$run.txt
         cp results_$run.txt ../
      else
         echo 'could not find hsoft or hfact file'
         exit
      endif
   endif
   cat results_$run.txt
   echo 'done'
endif
