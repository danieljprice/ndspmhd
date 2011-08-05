#!/usr/bin/perl
#
# swaps column order in my .ev files to be compatible with Matthew's grab output
#
use strict;
my $zero = 0.0;
while (<STDIN>) {
       my( $time, $ekin, $etherm, $emag, $epot, $etot, $linmom, undef ) = split ' ',$_;
       print "$time $time $time $etot $ekin $epot $etherm $zero $linmom \n";   
    }
