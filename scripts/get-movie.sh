#!/bin/bash
n=53;
rm -f splash_0*.png;
#splash -p prog -y 2 -x 1 -dev w.png `head -$n prog.filenames`;
#splash -p prog -y 3 -x 1 -dev gradw.png `head -$n prog.filenames`;
#splash -p prog -y 4 -x 1 -dev grgrw.png `head -$n prog.filenames`;
#nsplash -y 2 -x 1 -dev vx.png `head -$n splash.filenames`;
#nsplash -y 19 -x 1 -dev vx.png -p diff `head -$n splash.filenames`;
for i in `seq -w 0 $(( n - 1 ))`; do
    a="vx_00${i}.png";
    b=${a/vx/w};
    c=${a/vx/gradw};
    d=${a/vx/grgrw};
    out=${a/vx/splash};
    if [ -e $b ]; then
       pngtopnm $a > 1.ppm;
       pngtopnm $b > 2.ppm;
       pngtopnm $c > 3.ppm;
       pngtopnm $d > 4.ppm;
       pnmcat -lr 1.ppm 2.ppm > row1.ppm;
       pnmcat -lr 3.ppm 4.ppm > row2.ppm;
       echo "$a + $b + $c + $d -> $out";
       pnmcat -tb row1.ppm row2.ppm | pnmtopng > $out;
       echo "$a + $b + $c + $d -> $out";
       rm -f 1.ppm 2.ppm 3.ppm 4.ppm row1.ppm row2.ppm;
    fi
done
~/splash/scripts/movie.sh
