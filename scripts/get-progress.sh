#!/bin/bash
if [ -e progress ]; then
   echo "progress->progress.old";
   mv progress progress.old
fi
rm prog.filenames;
rm splash.filenames;
for x in gen? gen?? gen???; do
    num=${x/gen/};
    echo "$num `head -1 $x/ranked.list`" >> progress;
    best=`head -1 $x/ranked.list | cut -d' ' -f 1`;
    echo $x/$best.dat >> prog.filenames;
    echo "$x/$best/diff_00001.dat" >> splash.filenames;
done
