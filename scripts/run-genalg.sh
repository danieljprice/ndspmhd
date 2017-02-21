#!/bin/bash
#
# script to initiate and run kernel-breeding genetic algorithm sequence
#
# The plugin is a `run-child.sh` script that runs the child
# This script is specific to each calculation
##----------------
##!/bin/bash
# example run-child script
# relaxation
#     cp ../../unif.in .;
#     cp ../../rho.in .;
#     ../../2DSPMHD unif.in >& unif.log;
#     file=rho_00000.dat;
#     cp unif_00001.dat $file;
#     ../../2DSPMHD $file >& rho.log;
#     line=`nsplash -y 7 -x 1 -dev /null $file -p ../../splash | grep L2`;
#
##!/bin/bash
# example run-child script
# sod shock
#     cp ../../sod.in .;
#     cp ../../sod.shk .;
#     ../../1DSPMHD sod.in >& /dev/null;
#     file=sod_00001.dat;
#     line=`nsplash -y 6 -x 1 -dev /null $file -p ../../sod | grep L2`;
##!/bin/bash
# example run-child script
# 2nd deriv - diffusion problem
#     cp ../../diff.in .;
#     ../../3DSPMHD diff.in >& /dev/null;
#     rm diff_00000.dat diff.in kernel??.dat version diff.ev;
#     file=diff_00001.dat;
#     line=`nsplash -y 19 -x 1 -dev /null $file -p ../../diff | grep L2`;
##----------------
NDSPMHD_DIR=~/ndspmhd;
nkids=20;
#
# run a complete generation of children
#
run_gen() 
{
 rm -f score.list;
 for child in `seq -w 1 $nkids`; do
     if [ ! -d $child ]; then
        mkdir $child;
     fi
     #echo $child;
     cd $child;
     rm -f kernel100.dat;
     ln -s ../$child.dat kernel100.dat;
     #line=`../../run-child.sh`;
     #score=`echo $line | cut -d'=' -f 3 | cut -d'L' -f 1`;
     score=`../../run-child.sh`;
     mutfile="../$child.dat.mut";
     mut='';
     if [ -e $mutfile ]; then
        mut=`cat $mutfile`;
     fi
     echo $child $score $mut > score;
     cat score;
     echo $child $score $mut >> ../score.list;
     cd ..;
 done
 ../sort score.list > ranked.list;
}

#
# initiate the first generation
#
make_starting_gen()
{
  #for x in `seq -w 1 20`; do
  #   cp startlib/$x.dat .
  #   echo $x;
  #  ./breed cubic.dat cubic.dat $x.dat;
  #done
  cd startlib;
  run_gen;
  cd ..;
  ./breed startlib
}

#
# create a new generation from breeding the previous one
#
make_new_gen()
{
  dir=$1;
  cat $dir/ranked.list;
  ./breed $dir;
}

#
# preliminaries
#
if [ ! -e 01.dat ]; then
   make_starting_gen;
fi
if [ ! -e breed ]; then
   echo "compiling breed...";
   cd $NDSPMHD_DIR/utils; make breed; cd -; cp $NDSPMHD_DIR/utils/breed .;
fi
if [ ! -e sort ]; then
   echo "compiling sort...";
   cd $NDSPMHD_DIR/utils; make sort; cd -; cp $NDSPMHD_DIR/utils/sort .;
fi
#
# perform the loop over generations
#
for n in `seq 1 999`; do
    dir=gen$n;
    echo $dir;
    mkdir $dir;
    cd $dir;
    mv ../??.dat* .;
    run_gen;
    cd ..
    make_new_gen $dir;
done
