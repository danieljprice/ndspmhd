#!/usr/local/bin/perl
#
#
#
use strict;

my $n;
my $nruns;
my $nsteps;
my $startnum;
my $zero = 0;
my $nstart = 1;	## number to start from (either zero or 1)
my $ext = '.dat';

if ($#ARGV !=2) {
   die "Usage: $0 rootname nruns nstepsmax\n";
}

my ($rootname,$nruns,$nsteps) = @ARGV;

# call the multirun program to generate the input files
system "multi/multirun $rootname $nruns";

# write appropriate runnames to 'runname' and execute program
for ($n = $nstart;$n<=$nruns;$n++) {
    print "doing run $rootname$n \n";
    system "1DSPMHD $rootname$n";
}

# call the sametime program to write simultaneous data steps
system "utils/sametime $rootname $nruns $nsteps";

print "finished $rootname";

exit;
