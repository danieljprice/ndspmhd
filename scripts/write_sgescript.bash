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
echo '#$ -S /usr/bin/tcsh'                       ## run in tcsh
echo '#$ -cwd'					 ## run in current dir
echo '#$ -M dprice@ast.cam.ac.uk'		 ## mail me 
echo '#$ -m ae'					 ## (on abort or end)
echo '#$ -o '$PWD'/'$filename'.output'		 ## output from job in filename
echo '#$ -e '$PWD'/'$filename'.errors'		 ## errors from job in filename
echo '#$ -l s_vmem=30M'				 ## max memory requirements
echo '#$ -l h_vmem=30M'				 ##  "    "        "
echo '#$ -N' `/home/dprice/perl/randomword.pl`	 ## name job
echo $PWD'/1DSPMHD '$filename			 ## run the job
fi




