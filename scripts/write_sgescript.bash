#!/usr/bin/bash
#
# @(#) writes script for job submission to sun grid engine -- daniel price 27/10/03
#
if [ $# != 1 ]; then
   echo "Usage: $0 scriptname";
else
   filename=$1

#--------------------- print the script -------------------------------

echo '## Sun Grid Engine Script, created '`date`
echo '#$ -S /usr/bin/tcsh'
echo '#$ -cwd'
echo '#$ -M dprice@ast.cam.ac.uk'
echo '#$ -m ae'
echo '#$ -o '$PWD'/'$filename'.output'
echo '#$ -e '$PWD'/'$filename'.errors'
echo '#$ -l s_vmem=30M'
echo '#$ -l h_vmem=30M'
echo '#$ -N' `/home/dprice/perl/randomword.pl`
echo $PWD'/1DSPMHD '$filename
fi




