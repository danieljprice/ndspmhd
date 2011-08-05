#!/bin/tcsh
#
# @(#) makes new directory, copies input file and executable into it
# @(#)  and runs the job -- Daniel Price 27/10/03
#
if $# != 1 then
   echo Usage: $0 runname
else
   set runname=$1
   set dirname = ./$runname
   echo 'making new directory' $dirname
   mkdir $dirname
   echo 'copying files...'
   cp $runname.in $dirname
   cp 3DSPMHD $dirname
   cp defaults $dirname
   cd $dirname
   ln -s ../supersphplot supersphplot
   ln -s ../evsupersph evsupersph
   echo 'writing sge script '$runname\.sge
   ../scripts/write_sgescript3D.bash $runname > $runname\.sge
##   echo 'submitting job to grid'
##   qsub $runname\.sge   
endif
