#!/bin/bash
#
# @(#) writes script for job submission to sun grid engine -- daniel price 17/2/03
#
if [ $# != 2 ]; then
   echo "Usage: $0 nD scriptname";
else
   filename=$2
   sgename=$filename.sge
   jobname=`/home/dprice/scripts/randomword.pl`

#--------------------- print the script -------------------------------

echo "Writing SGE script to file $sgename: jobname '$jobname'"
echo '## Sun Grid Engine Script, created '`date` > $sgename
echo '#$ -S /usr/bin/tcsh'   >> $sgename                    ## run in tcsh
echo '#$ -cwd'			>> $sgename		 ## run in current dir
echo '#$ -M dprice@ast.cam.ac.uk'	>> $sgename	 ## mail me 
echo '#$ -m ae'			>> $sgename		 ## (on abort or end)
echo '#$ -o '$PWD'/'$filename'.output'	>> $sgename	 ## output from job in filename
echo '#$ -e '$PWD'/'$filename'.errors'	>> $sgename	 ## errors from job in filename
echo '#$ -l s_vmem=30M'			>> $sgename	 ## max memory requirements
echo '#$ -l h_vmem=30M'		>> $sgename		 ##  "    "        "
echo '#$ -N' $jobname >> $sgename	 ## name job
echo $PWD'/'$1'SPMHD '$filename		>> $sgename	 ## run the job
fi




