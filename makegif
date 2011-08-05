#!/usr/bin/tcsh
#
# @(#) Copies pgplot files to new directory, renames them
# @(#)  and invokes gifmerge -- Daniel Price 22/7/03
#
if $# != 1 then
   echo Usage: $0 runname
else
   set newdir=./GIFS/$1
   echo making directory $newdir
   mkdir $newdir
   echo moving files...
   mv pgplot.gif* $newdir
   cd $newdir
   ../fixgifs		# bash script that renames the files
   gifmerge -5 pgplot.gif* > ../$1\.gif
   cd ..
   xanim $1\.gif &
   echo "that's all folks"
endif
