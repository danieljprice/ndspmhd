#!/usr/bin/tcsh
# makes copies of scripts in main directory
#
cd ..
ln -s ./scripts/run1D.tcsh run1D
ln -s ./scripts/run2D.tcsh run2D
ln -s ./scripts/run3D.tcsh run3D
ln -s ./scripts/splitdat.pl splitdat
ln -s ./scripts/makegif.tcsh makemovie
ln -s ./scripts/mergedat.pl mergedat
ln -s ./scripts/domulti.pl domulti
ln -s ./scripts/doparallel.pl doparallel1D 
ln -s ./scripts/doparallel2D.pl doparallel2D
