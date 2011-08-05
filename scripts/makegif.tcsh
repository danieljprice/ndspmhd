#!/usr/bin/tcsh
#
# @(#) Copies pgplot files to new directory, renames them
# @(#)  and invokes gifmerge -- Daniel Price 22/7/03
#
if $# != 1 then
   echo Usage: $0 runname
else
   set newdir=/home/dprice/ndspmhd/MOVIES/$1
   echo making directory $newdir
   mkdir $newdir
   echo moving files...
   mv pgplot.gif* $newdir
   cd $newdir
   /home/dprice/ndspmhd/scripts/fixgifs.bash		# bash script that renames the files
## make an animated gif   
   gifmerge -5 pgplot.gif* > ../$1\.gif
## also make a .fli animation
   ls pgplot.gif* > filelist
   fbm2fli -rx 850 -ry 680 -ox 20 filelist ../$1\.fli
   cd ..
   xanim $1\.fli &
   echo "that's all folks"
endif
