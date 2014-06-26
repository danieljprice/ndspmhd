##-------------------------------------------------------------------##
##     Makefile for ND-SUPERSPMHD      			             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2003	 	     ##
##-------------------------------------------------------------------##
#EDITOR = nedit

intro:
	@echo "specify number of dimensions, e.g. 'make 1D' or try 'make install'"

install: 
	@echo 'Compiling 1D code'
	cd src;  make 1D
	#cd src; make clean; make 2D
	#cd src; make clean; make 25D
	#cd src; make clean; make 3D

1D:
	cd src; make 1D
	
2D:
	cd src; make 2D

3D:
	cd src; make 3D

clean:
	cd src; rm *.o *.mod

cleanall: clean
	rm ?DSPMHD

edit1D:
	cd src; make edit1D

edit2D:
	cd src; make edit2D

edit3D:
	cd src; make edit3D

make: checkeditor
	cd src; $(EDITOR) Makefile &
	
checkeditor:
   ifndef EDITOR
	@echo 'ERROR: need EDITOR variable set for "make edit" (e.g. export EDITOR=emacs)'
   endif

version=${shell cat docs/version}
dist:
	@echo "building tarball of ndspmhd version ${version}"
	git archive public --format=tar.gz --prefix=ndspmhd/ > ndspmhd-v${version}.tar.gz
