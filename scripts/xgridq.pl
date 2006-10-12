#!/usr/bin/perl
#
# script to query xgrid output in a useful format
#
# written by Daniel Price 22/9/06
#
use strict;

my $xgridauth='-hostname cytosine.ex.ac.uk -auth Kerberos';
my @joblist= split(/,\s+/,`xgrid $xgridauth -job list `);
my $format=" %4u | %25s | %9s | %19s | %19s | %19s \n";
my $formath=" %4s | %25s | %9s | %19s | %19s | %19s \n";
printf $formath,'ID','Script name','Status','Date Submitted','Date Started','Date Stopped';

my $nrun = 0;
my $nfinished = 0;
my $npending = 0;
foreach (@joblist) {
   my ($jobid) = m/(\d+)/;  
   my $attributes = `xgrid $xgridauth -job attributes -id $jobid`;
   my ($status) = $attributes =~ m/jobStatus = (\w+);/;
   my ($name) = $attributes =~ m/name = "(\S+)";/;
   my ($submitted) = $attributes =~ m/dateSubmitted = (.+)\s\+.+;/;
   my ($started) = $attributes =~ m/dateStarted = (.+)\s\+.+;/;
   my ($stopped) = $attributes =~ m/dateStopped = (.+)\s\+.+;/;
   if ($status eq 'Running') {
      $nrun = $nrun + 1;
   } elsif ($status eq 'Finished') {
      $nfinished = $nfinished + 1;
#      my $runtime = $
   } elsif ($status eq 'Pending') {
      $npending = $npending + 1;
   }
#   print "$attributes \n";
   printf $format, $jobid, $name, $status, $submitted, $started, $stopped;
}
print " $nrun jobs running; $npending jobs pending; $nfinished jobs finished\n";
