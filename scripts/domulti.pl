#!/bin/env perl
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
my $SPMHD = 'SPMHD';

if ($#ARGV !=3) {
   die "Usage: $0 nD rootname nruns nstepsmax\n";
}

my ($ndim,$rootname,$nruns,$nsteps) = @ARGV;

# make a new directory for the run
print "making new directory $rootname";
system "mkdir $rootname";
# copy files to this directory
system "cp multirun.in ./$rootname";
system "cp $ndim$SPMHD ./$rootname";
system "cp defaults ./$rootname";
system "cd $rootname; ln -s ../supersphplot ./supersphplot";
system "cd $rootname; ln -s ../evsupersph ./evsupersph";

# call the multirun program to generate the input files
system "cd $rootname; ../multi/multirun $rootname $nruns";

# write appropriate runnames to 'runname' and execute program
for ($n = $nstart;$n<=$nruns;$n++) {
    print "doing run $rootname$n \n";
    system "cd $rootname; ./$ndim$SPMHD $rootname$n";
}

# call the sametime program to write simultaneous data steps
system "cd $rootname; ../utils/sametime $rootname $nruns $nsteps";

print "finished $rootname";

exit;
