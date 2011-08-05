#!/bin/bash
#
# @(#) build script for SPLASH releases
# @(#) -- Daniel Price
#
splashdir=~/splash
if [ $# -ne 2 ]; then
   echo Usage: $0 version description
   echo "(NB: run cvs -q tag 'v1_x_x--`date +%d/%m/%y`' and modify splash.f90 *first* before doing build)"
   nedit $splashdir/splash.f90 &
else
   date=`date '+%d/%m/%y'`
   vernum=$1
   echo 'date = '$date
   echo 'version = '$vernum
   echo 'comment = '$2
   builddir='splash'
   webdir=$HOME'/web/splash'
   echo 'creating build in directory '/tmp/$builddir
#
#--update release history
#
#--add comment line to version history file (text)
   cd $HOME/Documents/splash-releases
   echo $vernum':' > version-$vernum
   echo $2 >> version-$vernum
   echo ' ' >> version-$vernum
   cat version-0* version-1.?.* version-1.??.* > version_history
   cp version_history $splashdir/docs
   cp version_history $webdir/download/
#--add comment to latex version for documentation
   echo $vernum' & '$date' & '$2 '\\' > versiontex-$vernum.tex
   cat versiontex-0*.tex versiontex-1.?.*.tex versiontex-*.??.*.tex > version_history_tex.tex
   tail version_history_tex.tex
#--commit the new version history files to the cvs repository
   cp version_history_tex.tex $splashdir/docs/
   cp version_history $splashdir/docs/
#--version number to build into the documentation
   echo $vernum > $splashdir/docs/version
   cd $splashdir/docs; git commit -m 'version '$vernum version_history_tex.tex version_history version
   git stash;
   git svn dcommit;
   git stash pop;
#
#--get current cvs copy from the repository and rename directory
#
   cd /tmp
   rm -rf /tmp/splash
   svn export https://svn-vre.its.monash.edu.au/mathsci/splash
#   mv ./splash $builddir
#   rmdir ndspmhd
   cd $builddir
#
#  add various things to the build directory
#
#--version file
   echo splash, version $vernum, built on `date` > VERSION
#--changelog
   cd $splashdir
   ~/cvsutils/git2cl.pl > ChangeLog
   cp ChangeLog /tmp/$builddir
   cp /tmp/$builddir/ChangeLog $webdir/download
#--tar file of directory
   cd /tmp
   tarfile=$builddir-$vernum.tar
   tar cf $tarfile $builddir
   gzip $tarfile
   echo 'copying '$tarfile' to download directory'
   cd $webdir/download
   mv *.tar* ./archive
   #--get rid of stupid .DS_Store files
   rm .DS_Store ./archive/.DS_Store ../userguide/.DS_Store
   cd /tmp
   cp $tarfile.gz $webdir/download
   cd $webdir/; rm splash.tar.gz; ln -s download/$tarfile.gz splash.tar.gz
   cp $tarfile.gz ~/Documents/splash-releases/
   rm $webdir/download/INSTALL*
   cp /tmp/$builddir/INSTALL* $webdir/download
#
#--check that build compiles correctly
#
   echo 'checking build...'
   cd /tmp/$builddir
   make clean
   make DEBUG=yes all
#
#  build the docs for the website
#
   echo 'building documentation...'
   cd $splashdir/docs
   cvs update
#   cd /tmp/$builddir/docs
   pdflatex splash
   bibtex splash >> /tmp/latex.output
   pdflatex splash >> /tmp/latex.output
   pdflatex splash >> /tmp/latex.output
#   dvips splash -o splash.ps
#   cp splash.ps $webdir/userguide
#   ps2pdf13 splash.ps
#   gzip splash.ps
#   cp splash.ps.gz $webdir/userguide
   cp splash.pdf $webdir/userguide
   echo 'building html documentation...'
   nedit $webdir/rightmenu.html $webdir/newslatest.html $webdir/newsclip.html $webdir/download.html $webdir/news.html &
   latex2html -local_icons -image_type gif splash.tex
   rm -r $webdir/userguide/html;
   mv splash $webdir/userguide/html;
   rm *.aux *.blg *.dvi *.log *.toc;

#  add style to html pages
   cd $webdir/userguide/html
   mv splash.css crap.css
   cat crap.css $HOME/web/dan.css > splash.css
   for x in *.html; do
       echo 'adding style to '$x;
       cat $x | sed 's/\<\/BODY\>/\<\/div\>&/g' | sed 's/\<BODY \>/&\<div id="wrap"\>/g' | sed 's/\<BODY\>/&\<div id="wrap"\>/g' > crap.html;
       mv crap.html $x;
   done
   nedit $webdir/userguide/html/index.html &

#  fix permissions (just in case)
#
   cd $HOME/web;
   chmod -R a+r *
   
fi
