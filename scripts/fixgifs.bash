#!/usr/bin/bash
# @(#) renames pgplot filenames so they are listed in the correct order
# @(#) (replaces _ with _0 in single then double digit filenames)
# @(#) NB: only works for < 1000 files

for x in *_?;
do if test $x != '*_?'; then echo $x ${x/_/_0}; mv $x ${x/_/_0}; fi;
done;

for x in *_??;
do if test $x != '*_??'; then echo $x ${x/_/_0}; mv $x ${x/_/_0}; fi;
done;

##for x in *_???;
##do if test $x != '_???'; then echo $x ${x/_/_0}; mv $x ${x/_/_0}; fi;
##done;
