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
my $SPMHD = 'SPMHD';

if ($#ARGV <2) {
   die "Usage: $0 nD rootname nruns [nstartfrom]\n";
}

my ($ndim,$rootname,$nruns,$nstart) = @ARGV;
$nstart = 1 unless defined $nstart; ## number to start from (either zero or 1)

# make a new directory for the run
print "making new directory $rootname";
system "mkdir $rootname";
# copy files to this directory
system "cp multirun.in ./$rootname";
system "cp $ndim$SPMHD ./$rootname";
system "cp defaults ./$rootname";
#system "cd $rootname; ln -s ../supersphplot ./supersphplot";
#system "cd $rootname; ln -s ../evsupersph ./evsupersph";
system "cd $rootname; ../scripts/writemake.tcsh $ndim > Makefile";

# call the multirun program to generate the input files
system "cd $rootname; ../multi/multirun $rootname $nruns";

# write appropriate runnames to 'runname' and execute program
for ($n = $nstart;$n<=$nruns;$n++) {
    print "doing run $rootname$n \n";
    print "copying $rootname\_init.dat $rootname$n\_init.dat \n";
    ##system "cd $rootname; cp $rootname\_init.dat $rootname$n\_init.dat; ./$ndim$SPMHD $rootname$n\_init.dat > $rootname$n.output";
    ##system "cd $rootname; ./$ndim$SPMHD $rootname$n\_init.dat > $rootname$n.output";
    system "cd $rootname; ./$ndim$SPMHD $rootname$n > $rootname$n.output < crap.in";
}

# call the sametime program to write simultaneous data steps
#system "cd $rootname; ../utils/sametime $rootname $nruns $nsteps";

print "finished $rootname \n";

exit;
