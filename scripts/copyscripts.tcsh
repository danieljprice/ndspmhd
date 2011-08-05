#!/bin/tcsh
# makes copies of scripts in main directory
#
cd ..
rm run splitdat makemovie mergedat domulti doparallel write_sgescript
ln -s ./scripts/mvdat.sh mvdat
ln -s ./scripts/run.tcsh run
ln -s ./scripts/splitdat.pl splitdat
ln -s ./scripts/makegif.tcsh makemovie
ln -s ./scripts/mergedat.pl mergedat
ln -s ./scripts/domulti.pl domulti
ln -s ./scripts/doparallel.pl doparallel
ln -s ./scripts/write_sgescript.bash write_sgescript
