#!/usr/bin/tcsh
# makes copies of scripts in main directory
#
set dirname='..'
ln -s run1D.tcsh $dirname/run1D
ln -s run2D.tcsh $dirname/run2D
ln -s run3D.tcsh $dirname/run3D
ln -s splitdat.pl $dirname/splitdat
ln -s makegif.tcsh $dirname/makemovie
ln -s mergedat.pl $dirname/mergedat
ln -s domulti.pl $dirname/domulti
ln -s doparallel.pl $dirname/doparallel 
