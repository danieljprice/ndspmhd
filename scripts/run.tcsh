#!/bin/tcsh
#
# @(#) makes new directory, copies input file and executable into it
# @(#)  and runs the job -- Daniel Price 17/2/04
#
if $# != 2 then
   echo Usage: $0 nD runname
else
   set runname=$2
   set dirname = ./$runname
   echo 'making new directory' $dirname
   mkdir $dirname
   echo 'copying files...'
   cp $runname.in $dirname
   set runfile=$1'SPMHD'
   echo copying executable $runfile
   cp $1'SPMHD' $dirname
   cp defaults $dirname
   cd $dirname
   ln -s ../supersphplot supersphplot
   ln -s ../evsupersph evsupersph
   echo 'writing a Makefile'
   ../scripts/writemake.tcsh $1 > Makefile
   echo 'writing sge script '$runname\.sge
   ../scripts/write_sgescript.bash $1 $runname > $runname\.sge
##   echo 'submitting job to grid'
##   qsub $runname\.sge   
endif
