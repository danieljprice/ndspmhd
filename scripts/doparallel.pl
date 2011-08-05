#!/bin/env perl
#
# farms out jobs to a list of machines, finding available CPU
#
use strict;

my $n;
my $nruns;
my $nsteps;
my $startnum;
my $zero = 0;
my $ext = '.dat';
my $SPMHD = 'SPMHD';

if ($#ARGV !=2) {
   die "Usage: $0 nD rootname nruns \n";
}

my ($ndim,$rootname,$nruns) = @ARGV;

my (@machines) = `cat machinelist`;
my $pwd = `pwd`;
chomp($pwd);
if ( $#machines < 1 ) {
   die "ERROR: no machines specified in file machinelist \n";
} elsif ( $#machines < $nruns ) {
   print "WARNING: not enough machines given to run all jobs \n"
}


# make a new directory for the run
print "making new directory $rootname \n";
system "mkdir $rootname";
# copy files to this directory
system "cp multirun.in ./$rootname";
system "cp $ndim$SPMHD ./$rootname";
system "cp defaults ./$rootname";
##system "cd $rootname; ln -s ../supersphplot ./supersphplot";
##system "cd $rootname; ln -s ../evsupersph ./evsupersph";
system "cd $rootname; ../scripts/writemake.tcsh $ndim > Makefile";

# call the multirun program to generate the input files
system "cd $rootname; ../multi/multirun $rootname $nruns";


my $machine;
my $n = 1;
my $jobsrun;

# loop through all available machines looking for spare CPU
foreach $machine (@machines) {
    chomp($machine);
    print "----------------- \n trying $machine \n";
    # get load average for this machine using uptime command
    my $loadav = `ssh $machine uptime | cut -f5 -d':' | cut -f1 -d','`;
    chomp($loadav);
    print " load average last 1 minute = $loadav";
    if ( $loadav < 0.75 ) {
       print " ...OK \n";
       # run the job
       print "running $rootname$n on machine $machine at nice +19\n";
       system "ssh $machine 'cd $pwd/$rootname; nice +19 ./$ndim$SPMHD $rootname$n > $rootname$n.output &' ";
       $n = $n + 1;
       if ($n > $nruns) {
          print "===========================================\n";
          print "\n Hurrah! all jobs successfully submitted \n";
          $jobsrun = $n - 1;
          print " $jobsrun jobs run \n";
          exit;
       }   
    } else {
       print " ...too busy \n";
    }
}

$jobsrun = $n - 1;
print "=======================================================\n";
print "WARNING: not enough machines available to run all jobs \n";
print " $jobsrun jobs run \n";
exit;
