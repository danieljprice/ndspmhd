#!/usr/local/bin/perl
#
# @(#) perl script to merge different SPH data files into one file 
# @(#) so that different simulations can be compared
#
use warnings;
use strict;
use IO::File;

my $n;
my $nstep;
my $nfile;
my $nfiles;
my $filename;
my $outputfile;
my $ncol; my $nrow;
my @infile;

## read arguments from command line

if ($#ARGV <=1) {		## $# is highest index in array
   system "what $0"; print "\n";
   die "Usage: $0 nstepsmax filename1 filename2 ... filenameN\n";
}

my ($nstepsmax, @filenames) = @ARGV;

## work out number of files from input

$nfiles = $#filenames + 1;

print "nfiles = $nfiles, steps = $nstepsmax, files = @filenames \n";

## specify name of output file

($outputfile = $filenames[0]) =~ s/\.dat.*/_merged.dat/;

print "merged file = $outputfile \n";

## open output file

open(OUTFILE,"> $outputfile") || die("can't open $outputfile") ;

## open all the input files

for ($nfile = 0;$nfile < $nfiles;$nfile++) {
    $infile[$nfile] = new IO::File( '<' . $filenames[$nfile] )
    	or die "Couldn't open '$filenames[$nfile]': $!";
}

## loop over all the steps, read a timestep from each input file and output
## them sequentially into the output file

for ($nstep = 1;$nstep<=$nstepsmax;$nstep++) {
       
    ##print " timestep number $nstep \n";
       
    for ($nfile = 0;$nfile < $nfiles;$nfile++) {
         
	print "$filenames[$nfile] :";
        ## read header line from this timestep
        my $line = readline $infile[$nfile];
	##print "header line = $line \n";

	## determine number of rows, columns in timestep
        my( undef, undef, $nrow, undef, undef, undef, undef, $ncol ) = split ' ',$line;
        print " step $nstep rows = $nrow columns = $ncol \n";
	
	## write header line to merged file
	print OUTFILE "$line";
	
	## read data from this timestep and write to output file
        for ($n = 1;$n<=$nrow;$n++) {
	   $line = readline $infile[$nfile];
           ##print "$line";
	   print OUTFILE "$line";
	}
    }          
}   

## close all input files
foreach my $handle ( @infile ) {
    close $handle or die "Couldn't close handle: $!";
}

print "Finished: merged file = $outputfile \n";

## close the output file
close(OUTFILE);
