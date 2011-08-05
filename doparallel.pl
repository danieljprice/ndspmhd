#!/usr/local/bin/perl
#
# does multiple runs, but submits each one to the sun grid engine separately
#
use strict;

my $n;
my $nruns;
my $nsteps;
my $startnum;
my $zero = 0;
my $ext = '.dat';

if ($#ARGV !=2) {
   die "Usage: $0 rootname nruns nstepsmax\n";
}

my ($rootname,$nruns,$nsteps) = @ARGV;

# call the multirun program to generate the input files
system "multi/multirun $rootname $nruns";

# write appropriate runnames to 'runname' and execute program
for ($n = 1;$n<=$nruns;$n++) {
    print "doing run $rootname$n \n";
    system "write_sgescript $rootname$n";
    system "qsub $rootname$n\.sge";
}

# call the sametime program to write simultaneous data steps
###system "utils/sametime $rootname $nruns $nsteps";

# put the final timestep as the current plot file
###    $startnum = 1;
###    open(THISPLOT,"> thisplot") || die("can't open thisplot") ;
    
###    if ($nsteps < 10) {
###       print ">>> writing filename = $rootname$startnum$ext$zero$nsteps to 'thisplot' \n";    
###       print THISPLOT "$rootname$startnum$ext$zero$nsteps";
###    }   
###    else {
###       print ">>> writing filename = $rootname$startnum$ext$nsteps to 'thisplot' \n";
###       print THISPLOT "$rootname$startnum$ext$nsteps";    
###    }
###    close(THISPLOT);


exit;
