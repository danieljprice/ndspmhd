#!/usr/bin/perl
#
#
#
use strict;

my $n;
my $nruns;
my $nsteps;
my $startnum;
my $zero = 0;
my $ext = '.dat';
my $dir = `pwd`;
my $SPMHD = 'SPMHD';

if ($#ARGV <2) {
   die "Usage: $0 nD rootname nruns [nstartfrom]\n";
}

my ($ndim,$rootname,$nruns,$nstart) = @ARGV;
$nstart = 1 unless defined $nstart; ## number to start from (either zero or 1)
my $rundir = "$dir'/'$rootname";

# make a new directory for the run
print "making new directory $rootname";
system "mkdir $rootname";
# copy files to this directory
system "cp multirun.in ./$rootname";
system "cp $ndim$SPMHD ./$rootname";
system "cp defaults ./$rootname";
system "cd $rootname; ../scripts/writemake.tcsh $ndim > Makefile";

# call the multirun program to generate the input files
system "cd $rootname; ~/ndspmhd/multi/multirun $rootname $nruns";

# write appropriate runnames to 'runname' and execute program
for ($n = $nstart;$n<=$nruns;$n++) {
    print "doing run $rootname$n \n";
##    print "copying $rootname\_init.dat $rootname$n\_init.dat \n";
##    system "cd $rootname; cp $rootname\_init.dat $rootname$n\_init.dat; ./$ndim$SPMHD $rootname$n\_init.dat > $rootname$n.output";
    ##system "echo 'cd $dir ./$ndim$SPMHD $rootname$n ' > $rootname/$rootname$n.csh";
    print "cd $rootname; xgrid -hostname cytosine.ex.ac.uk -auth Kerberos -in \$PWD -job submit ./$ndim$SPMHD $rootname$n";
}

# call the sametime program to write simultaneous data steps
#system "cd $rootname; ../utils/sametime $rootname $nruns $nsteps";

print "finished $rootname \n";

exit;
