#!/usr/local/bin/perl
#
#
#
use strict;

my $nplot;
my $nplots;
my $n;
my $dash;

if ($#ARGV !=1) {
   die "Usage: $0 rootname \n";
}
my ($rootname,$nplots) = @ARGV;

$dash = '_';

##system "rm supersphplot.log";

for ($nplot =1;$nplot<=$nplots;$nplot++) {
    print ">>> doing plot $nplot " ;
    system "cp defaults_$nplot defaults";
    system "supersphplot $rootname < supersphplot.txt ";
    system "mkdir ./GIFS/$rootname$dash$nplot";
    system "mv pgplot.gif* ./GIFS/$rootname$dash$nplot";
    system "cd ./GIFS/$rootname$dash$nplot";
# call bash script to rename the files
    system "../fixgifs";
# use gifmerge to create animation    
    system "gifmerge -5 pgplot.gif* > $rootname$dash$nplot.gif";
    system "cd ../../";
}

##foreach my $file ( glob '*_???' ) {
##	print "$file\n"
##}
##
##die;
