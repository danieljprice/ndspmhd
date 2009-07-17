#!/bin/tcsh
#
# @(#) writes a quick makefile for remaking the executable in the current 
# @(#) directory
#
if $# != 1 then
   echo Usage: $0 nD
else
   echo '#'
   echo '#--Makefile to remake the executable and copy to this directory'
   echo '#  Generated by '$0':' `date`
   echo '#'
   echo default':'$1
   echo
   echo '3D:'
   echo '	cd ~/ndspmhd/; make 3D; cd -; cp ~/ndspmhd/3DSPMHD .'
   echo '2D:'
   echo '	cd ~/ndspmhd/; make 2D; cd -; cp ~/ndspmhd/2DSPMHD .'
   echo '1D:'
   echo '	cd ~/ndspmhd/; make 1D; cd -; cp ~/ndspmhd/1DSPMHD .'
   echo '3DGR:'
   echo '	cd ~/ndspmhd/; make 3DGR; cd -; cp ~/ndspmhd/3DGRSPMHD .'
   echo '2DGR:'
   echo '	cd ~/ndspmhd/; make 2DGR; cd -; cp ~/ndspmhd/2DGRSPMHD .'
   echo '1DGR:'
   echo '	cd ~/ndspmhd/; make 1DGR; cd -; cp ~/ndspmhd/1DGRSPMHD .'
   echo 'plot:'
   echo '	cd ~/ndspmhd/plot; make ndspmhd'
   echo 'plotc:'
   echo '	cd ~/ndspmhd/plot; make clean'
   echo 'aplot:'
   echo '	cd ~/ndspmhd/plot; make ascii'
   echo 'edit:'
   echo '	cd ~/ndspmhd/; make edit'$1
   echo 'evplot:'
   echo '	cd ~/ndspmhd/evplot; make'
   echo 'evplotc:'
   echo '	cd ~/ndspmhd/evplot; make clean'
   echo 'moddump:'
   echo '	cd ~/ndspmhd/src; make moddump'$1
   echo 'moddumpc:'
   echo '	cd ~/ndspmhd/src; make clean; make moddump'$1
   echo 'make:'
   echo '	cd ~/ndspmhd/src; nedit Makefile &'
   echo 'clean:'
   echo '	cd ~/ndspmhd/; make clean'
endif
