#!/bin/bash
#
# @(#) writes script for job submission to sun grid engine -- daniel price 17/2/03
#
if [ $# != 2 ]; then
   echo "Usage: $0 nD scriptname";
else
   filename=$2

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
echo '#$ -N' `/home/dprice/scripts/randomword.pl`	 ## name job
echo $PWD'/'$1'SPMHD '$filename			 ## run the job
fi




