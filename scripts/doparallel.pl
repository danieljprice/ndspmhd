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
my $njobsrun = 1;
my $jobsrun;
my $ncpu = 1;
my $njobspermachine = 3;
my $n;

# loop through all available machines looking for spare CPU
foreach $machine (@machines) {
    chomp($machine);
    print "----------------- \n trying $machine \n";
    # get a list of CPU in use using ps and add it all up to get the total
    my @usages = `ssh $machine ps -eo %cpu | sed 1d `;
       # add up total cpu percentage currently in use
       my $cpuinuse = 0.0;
       for ($ncpu = 0;$ncpu<=$#usages;$ncpu++) {
           $cpuinuse = $cpuinuse + @usages[$ncpu];
       }
    print "cpu percentage in use = $cpuinuse \n";
#    get load average for this machine using uptime command
#    my $loadav = `ssh $machine uptime | cut -f5 -d':' | cut -f1 -d','`;
#    chomp($loadav);
#    print " load average last 1 minute = $loadav";
    if ( $cpuinuse < 50.0 ) {
       print " ...OK \n";
       # run the job
       for ($n = 1;$n<=$njobspermachine;$n++){
           print "running $rootname$njobsrun on machine $machine at nice +19\n";
           system "ssh $machine 'cd $pwd/$rootname; nice +19 ./$ndim$SPMHD $rootname$njobsrun > $rootname$njobsrun.output &' ";
           $njobsrun = $njobsrun + 1;
           if ($njobsrun > $nruns) {
              print "===========================================\n";
              print "\n Hurrah! all jobs successfully submitted \n";
              $jobsrun = $njobsrun - 1;
              print " $jobsrun jobs run \n";
              exit;
           }  
       } 
    } else {
       print " ...too busy \n";
    }
}

$jobsrun = $njobsrun - 1;
print "=======================================================\n";
print "WARNING: not enough machines available to run all jobs \n";
print " $jobsrun jobs run \n";
exit;
