#!/bin/env perl
#
# parses COMPILER file for compiler flags and changes them
# appropriately in all of the Makefiles
#
# There are probably much better (shorter) ways of doing this, but it works!
#
# Flags are only changed if they are UNCOMMENTED in the Makefiles
#
use strict;
my $flag;
my $flagsetting;
my $f90C;
my $f90flags;
my $ldflags;
my $compilerfile = 'COMPILER';
my $line;
my $crap;
my @makefiles = ("./utils/Makefile","./evplot/Makefile","./plot/Makefile","./multi/Makefile","./src/Makefile");
#
# open the file COMPILER and read the flag settings
#
open(FILE,"< $compilerfile") || die("can't open file $compilerfile") ;

   foreach $line (<FILE>) {
       if ($line!~/^#/) {
          my( $flag, $flagsetting ) = split '=',$line;

          if ($flag=~/^(F90C)/) {
	     $f90C = $flagsetting;
	     #print "F90C = $f90C";
          } elsif ($flag=~/^(F90FLAGS)/) {
	     $f90flags = $flagsetting;
	     #print "F90FLAGS = $f90flags"
          } elsif ($flag=~/^(LDFLAGS)/) {
	     $ldflags = $flagsetting;
	     #print "LDFLAGS = $ldflags"
	  }       
       }
   }   
close(FILE);
#
# check that all appropriate flags have been input
#
print "------------- new flag settings ------------------- \n";
print "F90C = $f90C";
print "F90FLAGS = $f90flags";
print "LDFLAGS = $ldflags";
#
# Now go through all the Makefiles and replace the appropriate flags
#
my $makefile;
my $date = `date`;
print $date;

foreach $makefile(@makefiles) {

   print "------------------------------------------------------ \n";   
   ## make backup copy of each file
   my $oldmakefile = "$makefile.old";
   print "copying $makefile to $oldmakefile \n";
   system "cp $makefile $oldmakefile" || die("can't write $oldmakefile");

   open(OLDMAKEFILE,"< $oldmakefile") || die("can't open $oldmakefile") ;
   open(MAKEFILE,"> $makefile") || die("can't open $makefile") ;

   print "changing flag settings in $makefile...";
      foreach $line (<OLDMAKEFILE>) {
          if ($line!~/^#/) {
             my( $flag, $flagsetting ) = split '=',$line;

             if ($flag=~/^(F90C)/) {
	        print MAKEFILE "F90C = $f90C";
             } elsif ($flag=~/^(F90FLAGS)/) {
	        print MAKEFILE "F90FLAGS = $f90flags";
             } elsif ($flag=~/^(LDFLAGS)/) {
	        print MAKEFILE "LDFLAGS = $ldflags";
	     } else {
	        print MAKEFILE $line;
	     }       
          } else {
	     print MAKEFILE $line;
	  }
      }   
      print " done\n";
   close(OLDMAKEFILE);
   close(NEWMAKEFILE);

}
