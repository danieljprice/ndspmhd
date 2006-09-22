#!/usr/bin/perl
# adds up the L2 errors from supersphplot output; calculates average
use strict;
use warnings;

use List::Util qw(sum);
my $file = ' ';
my $avg = 0;
my $avg1 = 0;
my $avg2 = 0;
if ($#ARGV < 0) {
   print "script which parses supersphplot output for L2 errors \n";
   print "and calculates the average. Written by D. Price. \n\n";
   die "Usage: $0 filename(s) \n";
}

foreach $file (@ARGV) {
   open my $fh, '<', $file or die "Can't open $file: $!";
   my @errors;
   my @errors1;
   my @errors2;

   while ( <$fh> ) {
      if ( my ($error) = m/L2 error\s+=\s+(\S*)\s+/ ) {
         #print "error = $error \n";
         push @errors, $error;
      }
      if ( my ($error1) = m/L2 err \(1\)\s+=\s+(\S*)\s+/ ) {
         #print "error = $error \n";
         push @errors1, $error1;
      }
      if ( my ($error2) = m/L2 err \(2\)\s+=\s+(\S*)\s+/ ) {
         #print "error = $error \n";
         push @errors2, $error2;
      }
   }

   my $nerrors = scalar(@errors);
   my $nerrors1 = scalar(@errors1);
   my $nerrors2 = scalar(@errors2);
   if ($nerrors > 0) {
      $avg = sum(@errors) / $nerrors;
   }
   else {
      $avg = 0;   
   }
   if ($nerrors1 > 0) {
      $avg1 = sum(@errors1) / $nerrors1;
   }
   else {
      $avg1 = 0;   
   }
   if ($nerrors2 > 0) {
      $avg2 = sum(@errors2) / $nerrors2;
   }
   else {
      $avg2 = 0;   
   }

   print "$file: Average of $nerrors errors: $avg : $avg1 : $avg2 : (nerrors = $nerrors1, $nerrors2 )\n";
}

#print "$_\n" for @errors;
