#!/bin/env perl
#
# does multiple runs, but submits each one to the sun grid engine separately
#
use strict;

my $n;
my $n1;
my $nruns;
my $nsteps;
my $startnum;
my $zero = 0;
my $scriptfile;
my $slash = "/";
my $cpucheck = 'checkcpu.pl';

if ($#ARGV !=1) {
   die "Usage: $0 rootname nruns \n";
}

my ($rootname,$nruns) = @ARGV;

# save the location of this directory
my $thisdir = `pwd`;
$cpucheck = "$thisdir$slash$cpucheck";

# make a new directory for the run
print "making new directory $rootname \n";
system "mkdir $rootname";
my $newdir = `cd $rootname; pwd`;
print "$newdir ";
# copy files to this directory
system "cp multirun.in ./$rootname";
system "cp 1DSPMHD ./$rootname";
system "cp defaults ./$rootname";
system "ln -s ./plot/supersphplot ./$rootname/supersphplot";
system "ln -s ./evsupersph ./$rootname/evsupersph";

# call the multirun program to generate the input files
system "cd $rootname; ../multi/multirun $rootname $nruns";

# read in a list of machines to run jobs on
my @machines = `cat listofmachines`;
if ($#machines < $nruns) {die 'not enough machines to run on';}

# write appropriate runnames to 'runname' and execute program
for ($n = 1;$n<=$nruns;$n++) {
    $n1 = $n - 1;
# run the job remotely if enough cpu space available;
    my $machine = @machines[$n1];
    print "$machine why the line break";
    my $cpuinuse = `ssh cornfed /home/dprice/ndspmhd/checkcpu.pl`;
    print "cpu percentage in use = $cpuinuse";
    if ( $cpuinuse < 50.0 ) {print "plenty of space available..."};
    ##system "ssh $machine nice 19 ./1DSPMHD $rootname$n \>\>& $rootname$n.output \& \n";
# run this script
    print "doing run $rootname$n on machine @machines[$n1]\n";
##    system "./$rootname/$scriptfile"    
}

# call the sametime program to write simultaneous data steps
###system "utils/sametime $rootname $nruns $nsteps";

exit;
