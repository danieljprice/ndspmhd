#!/bin/tcsh
#
# @(#) build script for SPLASH releases
# @(#) -- Daniel Price
#
if $# != 2 then
   echo Usage: $0 version description
   echo "(NB: run cvs -q tag 'v1_x_x--`date +%d/%m/%y`' and modify splash.f90 *first* before doing build)"
   nedit ~/ndspmhd/plot/splash.f90 &
else
#
#--get current cvs copy from the repository and rename directory
#
   set date=`date '+%d/%m/%y'`
   set vernum=$1
   echo 'date = '$date
   echo 'version = '$vernum
   echo 'comment = '$argv[2]
   set builddir='splash'
   echo 'creating build in directory '/tmp/$builddir
   cd /tmp
   rm -r /tmp/ndspmhd
   rm -r /tmp/splash
   cvs export -Dtoday ndspmhd/plot
   mv ./ndspmhd/plot $builddir
   rmdir ndspmhd
   cd $builddir
#
#  add various things to the build directory
#
#--version file
   echo splash, version $vernum, built on `date` > VERSION
#--version number to build into the documentation
   echo $vernum > ./docs/version
#--changelog
   cd ~/ndspmhd/plot
   cvs2cl -t
   cp ChangeLog /tmp/$builddir
   cp /tmp/$builddir/ChangeLog ~/web/splash/download
#--tar file of directory
   cd /tmp
   set tarfile = $builddir-$vernum.tar
   tar cf $tarfile $builddir
   gzip $tarfile
   echo 'copying '$tarfile' to download directory'
   cd ~/web/splash/download
   mv *.tar* ./archive
   #--get rid of stupid .DS_Store files
   rm .DS_Store ./archive/.DS_Store ../userguide/.DS_Store
   cd /tmp
   cp $tarfile.gz ~/web/splash/download
   cd ~/web/splash/; rm splash.tar.gz; ln -s download/$tarfile.gz splash.tar.gz
   cp $tarfile.gz ~/Documents/splash-releases/
   cp $builddir/INSTALL* ~/web/splash/download
#--add comment line to version history file (text)
   cd ~/Documents/splash-releases
   echo $vernum':' > version-$vernum
   echo $argv[2] >> version-$vernum
   echo ' ' >> version-$vernum
   cat version-* > version_history
   cp version_history /tmp/$builddir/docs
   cp version_history ~/web/splash/download/
#--add comment to latex version for documentation
   echo $vernum" & "$date" & "$argv[2] "\\" > versiontex-$vernum.tex
   cat versiontex-*.tex > version_history_tex.tex
   cp version_history_tex.tex /tmp/$builddir/docs
#
#--check that build compiles correctly
#
   echo 'checking build...'
   cd /tmp/$builddir
   make clean
   make all
#
#  build the docs for the website
#
   echo 'building documentation...'
   cd /tmp/$builddir/docs
   latex splash
   bibtex splash >> /tmp/latex.output
   latex splash >> /tmp/latex.output
   latex splash >> /tmp/latex.output
   dvips splash -o splash.ps
   cp splash.ps ~/web/splash/userguide
   ps2pdf13 splash.ps
   gzip splash.ps
   cp splash.ps.gz ~/web/splash/userguide
   cp splash.pdf ~/web/splash/userguide
   echo 'building html documentation...'
   nedit ~/web/splash/index.html &
   latex2html splash.tex
   rm -r ~/web/splash/userguide/html
   mv splash ~/web/splash/userguide/html
   rm *.aux *.blg *.dvi *.log *.toc
   cd ~/web
   chmod -R a+r *
   
endif
