#!/bin/env perl
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

if ($#ARGV !=1) {
   die "Usage: $0 rootname nruns \n";
}

my ($rootname,$nruns) = @ARGV;

# make a new directory for the run
print "making new directory $rootname";
system "mkdir $rootname";
# copy files to this directory
system "cp multirun.in ./$rootname";
system "cp 1DSPMHD ./$rootname";
system "cp defaults ./$rootname";
system "ln -s ./plot/supersphplot ./$rootname/supersphplot";
system "ln -s ./evsupersph ./$rootname/evsupersph";

# call the multirun program to generate the input files
system "cd $rootname; ../multi/multirun $rootname $nruns";

# write appropriate runnames to 'runname' and execute program
for ($n = 1;$n<=$nruns;$n++) {
    print "doing run $rootname$n \n";
    system "cd $rootname; ../scripts/write_sgescript1D.bash $rootname$n > ./$rootname$n\.sge";
    system "cd $rootname; qsub $rootname$n\.sge";
}

# call the sametime program to write simultaneous data steps
###system "utils/sametime $rootname $nruns $nsteps";

exit;
