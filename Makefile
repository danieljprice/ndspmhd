##-------------------------------------------------------------------##
##     Makefile for ND-SUPERSPMHD      			             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2003	 	     ##
##-------------------------------------------------------------------##

SAVEDIR = ../version_control/`cat version`_`date +%d_%m_%Y`
EDITOR = nedit

install: 
	echo 'Making copies of useful scripts...'
	cd scripts; ./copyscripts.tcsh
	echo 'Compiling 1D code'
	cd src;  make 1D
	#cd src; make clean; make 2D
	#cd src; make clean; make 25D
	#cd src; make clean; make 3D
all:	
	echo 'Making copies of useful scripts...'
	cd scripts; ./copyscripts.tcsh
	echo 'Compiling plotting utility...'
	cd plot; make dansph; make mrbsph	
	echo 'Compiling .ev plotting utility...'
	cd evplot; make
	echo 'Compiling other utilities...'  
	cd utils; make sametime
	cd multi; make initial; make multirun
	echo 'Compiling 1D code'
	cd src;  make 1D
	
plotsph:
	cd plot; make

plotev:
	cd evplot; make
	
multirun:
	cd multi; make

1D:
	cd src; make 1D
	
2D:
	cd src; make 2D

25D:
	cd src; make 25D

3D:
	cd src; make 3D

clean:
	cd src; rm *.o *.mod

edit1D:
	cd src; make edit1D

edit2D:
	cd src; make edit2D
	
edit25D:
	cd src; make edit25D

edit3D:
	cd src; make edit3D

editev:
	cd evplot; $(EDITOR) evsupersph.f &

scripts:	## makes symbolic links to scripts
	ln -s run2D ./scripts/run2D.tcsh
	ln -s run3D ./scripts/run3D.tcsh
make:
	cd src; $(EDITOR) Makefile &

tag:
	cvs -q tag `cat version`_`date +%d_%m_%Y`
	
save:
	mkdir $(SAVEDIR)
	mkdir $(SAVEDIR)/src
	cd src; make targz
	mv ./src/*.tar.gz $(SAVEDIR)/src
	cp *.in $(SAVEDIR)
	cp *.pl $(SAVEDIR)
	mkdir $(SAVEDIR)/multi
	cp ./multi/*.f90 $(SAVEDIR)/multi
	cp ./multi/Makefile $(SAVEDIR)/multi
	cd plot; make targz
	mkdir $(SAVEDIR)/plot
	mv ./plot/*.tar.gz $(SAVEDIR)/plot
	mkdir $(SAVEDIR)/utils
	cp ./utils/*.f90 $(SAVEDIR)/utils
	cp -r evplot $(SAVEDIR)
	cp unpackall $(SAVEDIR)
