#!/bin/tcsh
#
# @(#) makes new directory, copies input file and executable into it
# @(#)  and runs the job -- Daniel Price 17/2/04
#
if $# != 2 then
   echo Usage: $0 version description
   echo "(NB: run cvs -q tag 'v1_x_x--`date +%d/%m/%y`' and modify supersphplot.f90 *first* before doing build)"
   nedit ~/ndspmhd/plot/supersphplot.f90 &
else
#
#--get current cvs copy from the repository and rename directory
#
   set date=`date '+%d/%m/%y'`
   set vernum=$1
   echo 'date = '$date
   echo 'version = '$vernum
   echo 'comment = '$argv[2]
   set builddir='supersphplot'
   echo 'creating build in directory '/tmp/$builddir
   cd /tmp
   rm -r /tmp/ndspmhd
   rm -r /tmp/supersphplot
   cvs export -Dtoday ndspmhd/plot
   mv ./ndspmhd/plot $builddir
   rmdir ndspmhd
   cd $builddir
#
#  add various things to the build directory
#
#--version file
   echo SUPERSPHPLOT, version $vernum, built on `date` > VERSION
#--version number to build into the documentation
   echo $vernum > ./docs/version
#--changelog
   cd ~/ndspmhd/plot
   cvs2cl -t
   cp ChangeLog /tmp/$builddir
   cp /tmp/$builddir/ChangeLog ~/web/supersphplot/download
#--tar file of directory
   cd /tmp
   set tarfile = $builddir-$vernum.tar
   tar cf $tarfile $builddir
   gzip $tarfile
   echo 'copying '$tarfile' to download directory'
   cd ~/web/supersphplot/download
   mv *.tar* ./archive
   #--get rid of stupid .DS_Store files
   rm .DS_Store ./archive/.DS_Store ../userguide/.DS_Store
   cd /tmp
   cp $tarfile.gz ~/web/supersphplot/download
   cp $tarfile.gz ~/supersphplot-releases/
   cp $builddir/INSTALL* ~/web/supersphplot/download
#--add comment line to version history file (text)
   cd ~/supersphplot-releases
   echo $vernum':' > version-$vernum
   echo $argv[2] >> version-$vernum
   echo ' ' >> version-$vernum
   cat version-* > version_history
   cp version_history /tmp/$builddir/docs
#--add comment to latex version for documentation
   echo $vernum" & "$date" & "$argv[2] "\\" > versiontex-$vernum.tex
   cat versiontex-*.tex > version_history_tex.tex
   cp version_history_tex.tex /tmp/$builddir/docs
#
#  build the docs for the website
#
   echo 'building documentation...'
   cd /tmp/$builddir/docs
   latex supersphplot
   bibtex supersphplot >> /tmp/latex.output
   latex supersphplot >> /tmp/latex.output
   latex supersphplot >> /tmp/latex.output
   dvips supersphplot -o supersphplot.ps
   cp supersphplot.ps ~/web/supersphplot/userguide
   ps2pdf13 supersphplot.ps
   gzip supersphplot.ps
   cp supersphplot.ps.gz ~/web/supersphplot/userguide
   cp supersphplot.pdf ~/web/supersphplot/userguide
   echo 'building html documentation...'
   nedit ~/web/supersphplot/index.html &
   latex2html supersphplot.tex
   rm -r ~/web/supersphplot/userguide/html
   mv supersphplot ~/web/supersphplot/userguide/html
   rm *.aux *.blg *.dvi *.log *.toc
   cd ~/web
   chmod -R a+r *
   
endif
