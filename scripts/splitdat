#!/usr/local/bin/perl
#
# @(#) perl script to split SPH data file into separate files containing 
# @(#) one timestep each (plus header line)
#
use strict;

my $n;
my $nstepsmax;
my $nstep;
my $filename;
my $ncol; my $nrow;

if ($#ARGV !=1) {
   die "Usage: $0 filename nstepsmax\n";
}

my ($filename,$nstepsmax) = @ARGV;


open(FILE,"< $filename") || die("can't open $filename") ;

   for ($nstep = 1;$nstep<=$nstepsmax;$nstep++) {
       my $line = <FILE>;
       my( undef, undef, $nrow, undef, undef, undef, undef, $ncol ) = split ' ',$line;
       print "step $nstep rows = $nrow columns = $ncol \n";
   
       open(STEPFILE,"> $filename\_$nstep") || die(" can't open $filename\_$nstep");
   
       print "writing step $nstep to file $filename\_$nstep\n";
       print "$line";
       print STEPFILE "$line";
   
       for ($n = 1;$n<=$nrow;$n++) {
           my $line = <FILE>;
           ###print "$line";
	   print STEPFILE "$line";
       }
  
       close(STEPFILE);
       
   }   
close(FILE);
