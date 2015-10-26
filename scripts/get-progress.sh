#!/bin/bash
if [ -e progress ]; then
   echo "progress->progress.old";
   mv progress progress.old
fi
if [ $# -ge 0 ]; then
   prefix=$1;
   echo "filename is ${prefix}_00001.dat";
else
   prefix='diff';
fi
ntop=1;
rm prog.filenames;
rm splash.filenames;
for x in gen? gen?? gen???; do
    num=${x/gen/};
    echo "$num `head -1 $x/ranked.list`" >> progress;
    if [ $ntop > 1 ]; then
       best=`head -1 $x/ranked.list | cut -d' ' -f 1`;
       echo $x/$best.dat >> prog.filenames;
#      echo "$x/$best/${prefix}_00001.dat" >> splash.filenames;
       echo $x/$best/errors-cp.out >> splash.filenames;    
    else
       best=`head -1 $x/ranked.list | cut -d' ' -f 1`;
       echo $x/$best.dat >> prog.filenames;
#      echo "$x/$best/${prefix}_00001.dat" >> splash.filenames;
       echo $x/$best/errors-cp.out >> splash.filenames;
    fi
done
