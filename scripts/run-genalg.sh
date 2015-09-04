run_gen() 
{
 rm -f score.list;
 for child in `seq -w 1 20`; do
     if [ ! -d $child ]; then
        mkdir $child;
     fi
     #echo $child;
     cd $child;
     rm -f kernel100.dat;
     ln -s ../$child.dat kernel100.dat;
# relaxation
#     cp ../../unif.in .;
#     cp ../../rho.in .;
#     ../../2DSPMHD unif.in >& unif.log;
#     file=rho_00000.dat;
#     cp unif_00001.dat $file;
#     ../../2DSPMHD $file >& rho.log;
#     line=`nsplash -y 7 -x 1 -dev /null $file -p ../../splash | grep L2`;
# sod shock
#     cp ../../sod.in .;
#     cp ../../sod.shk .;
#     ../../1DSPMHD sod.in >& /dev/null;
#     file=sod_00001.dat;
#     line=`nsplash -y 6 -x 1 -dev /null $file -p ../../sod | grep L2`;
# 2nd deriv - diffusion problem
     cp ../../diff.in .;
     ../../3DSPMHD diff.in >& /dev/null;
     rm diff_00000.dat diff.in kernel??.dat version diff.ev;
     file=diff_00001.dat;
     line=`nsplash -y 19 -x 1 -dev /null $file -p ../../diff | grep L2`;
     errL2=`echo $line | cut -d'=' -f 3 | cut -d'L' -f 1`;
     errLinf=`echo $line | cut -d'=' -f 4 | cut -d'(' -f 1`;
     echo $child $errL2 $errLinf > score;
     cat score;
     echo $child $errL2 $errLinf >> ../score.list;
     cd ..;
 done
 ../sort score.list > ranked.list;
}

make_starting_gen()
{
  #for x in `seq -w 1 20`; do
  #   cp startlib/$x.dat .
  #   echo $x;
  #  ./breed cubic.dat cubic.dat $x.dat;
  #done
  ./breed startlib
}

make_new_gen()
{
  dir=$1;
  echo "doing nothing";
  cat $dir/ranked.list;
  ./breed $dir;
}

if [ ! -e 01.dat ]; then
   make_starting_gen;
fi
for n in `seq 1 999`; do
    dir=gen$n;
    echo $dir;
    mkdir $dir;
    cd $dir;
    mv ../??.dat .;
    run_gen;
    cd ..
    make_new_gen $dir;
done
