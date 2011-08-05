#!/usr/bin/env perl
#
# does multiple runs, but submits each one to the sun grid engine separately
#
use strict;

## get list of cpu percentages using ps
my @usages = `ps -eo %cpu | sed 1d `;
#print "@usages \n";

# add up total cpu percentage currently in use
my $sum = 0.0;
my $n;
for ($n = 0;$n<=$#usages;$n++) {
    $sum = $sum + @usages[$n];
}
print "$sum \n";

exit;
