#!/bin/bash
#
# utility to check the CVS root
#
if [ "$#" -ne 0 ]; then
   echo Warning: arguments ignored
fi

rootdir=$PWD;

for dir in `find . -type d`;
do 
   cd $dir;      
   if [ -d CVS ]; then
      cd CVS;
      echo $dir : CVS root is `cat Root`;
      cd .. ;
   fi
   cd $rootdir
done
