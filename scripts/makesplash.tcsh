#!/bin/tcsh
#
# @(#) builds SPLASH on all machines to which I have access
#
echo '------ obelix -------'
#ssh obelix.astro.ex.ac.uk "cd ndspmhd/plot; cvs update; make all"
echo '---- mac cluster ----'
ssh djp212@witch.ex.ac.uk "cd ndspmhd/plot; cvs update; make all"
echo '------ ukaff1a ------'
ssh -f ukaff1a.star.le.ac.uk "cd plot; cvs update; setenv PARALLEL 'yes'; make all"
