#!/bin/bash
# 
# utility to change the CVS root of an already checked-out code version
# changes the root in all subdirectories also
#
if [ "$#" -ne 1 ]; then
   echo Usage: $0 newcvsroot
else
   newcvsroot=$1;
   rootdir=$PWD;

   for dir in `find . -type d`;
   do 
      cd $dir;      
      if [ -d CVS ]; then
         cd CVS;
         echo $dir : old root was `cat Root`;
         echo $newcvsroot > Root ;
         cd .. ;
      ##else
      ##   echo no CVS directory found;
      fi
      cd $rootdir
   done
   echo NEW CVS ROOT IN ALL DIRS IS $newcvsroot
fi
