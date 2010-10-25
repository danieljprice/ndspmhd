##-------------------------------------------------------------------##
##     Makefile for ND-SUPERSPMHD      			             ##
##     Written by Daniel Price					     ##
##     Institute of Astronomy, Cambridge, UK, 2003	 	     ##
##-------------------------------------------------------------------##

SAVEDIR = ../version_control/`cat version`_`date +%d_%m_%Y`
#EDITOR = nedit

intro:
	@echo "specify number of dimensions, e.g. 'make 1D' or try 'make install'"

install: 
	@echo 'Compiling 1D code'
	cd src;  make 1D
	#cd src; make clean; make 2D
	#cd src; make clean; make 25D
	#cd src; make clean; make 3D
all:	
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

3D:
	cd src; make 3D

1DGR:
	cd src; make 1DGR
	
2DGR:
	cd src; make 2DGR

3DGR:
	cd src; make 3DGR

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

edit1DGR:
	cd src; make edit1D

edit2DGR:
	cd src; make edit2D

edit3DGR:
	cd src; make edit3D

editev: checkeditor
	cd evplot; $(EDITOR) evsupersph.f &

scripts:	## makes symbolic links to scripts
	ln -s run2D ./scripts/run2D.tcsh
	ln -s run3D ./scripts/run3D.tcsh
make: checkeditor
	cd src; $(EDITOR) Makefile &
	
makes: checkeditor
	$(EDITOR) ./src/Makefile ./evplot/Makefile ./plot/Makefile \
	./utils/Makefile ./multi/Makefile &

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
	mkdir $(SAVEDIR)/utils
	cp ./utils/*.f90 $(SAVEDIR)/utils
	cp -r evplot $(SAVEDIR)
	cp unpackall $(SAVEDIR)

checkeditor:
   ifndef EDITOR
	@echo 'ERROR: need EDITOR variable set for "make edit" (e.g. export EDITOR=emacs)'
   endif
