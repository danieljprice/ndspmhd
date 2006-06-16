#!/usr/bin/perl
#
#
#
use strict;

my $n;
my $nruns;
my $nsteps;
my $startnum;
my $zero = 0;
my $ext = '.dat';
my $dir = `pwd`;
my $SPMHD = 'SPMHD';
my $ndspmhddir = '~/ndspmhd';

if ($#ARGV <2) {
   die "Usage: $0 nD rootname nruns [nstartfrom]\n";
}

my ($ndim,$rootname,$nruns,$nstart) = @ARGV;
$nstart = 1 unless defined $nstart; ## number to start from (either zero or 1)
my $rundir = "$dir/$rootname";
print "rundir is $rundir \n";

# make a new directory for the run
print "making new directory $rootname \n";
system "mkdir $rootname" || die "error creating run directory";

# copy files to this directory
system "cp multirun.in ./$rootname" || die "ERROR: multirun.in not found/not readable";
system "cp $ndspmhddir/$ndim$SPMHD ./$rootname";
system "cd $rootname; $ndspmhddir/scripts/writemake.tcsh $ndim > Makefile";

# call the multirun program to generate the input files
system "cd $rootname; $ndspmhddir/multi/multirun $rootname $nruns";

# write appropriate runnames to 'runname' and execute program
for ($n = $nstart;$n<=$nruns;$n++) {
    print "doing run $rootname$n \n";
##    print "copying $rootname\_init.dat $rootname$n\_init.dat \n";
##    system "cd $rootname; cp $rootname\_init.dat $rootname$n\_init.dat; ./$ndim$SPMHD $rootname$n\_init.dat > $rootname$n.output";
    system "cd $rootname; rm $rootname$n.csh; echo '#!'/bin/tcsh > $rootname$n.csh; echo cd \$PWD >> $rootname$n.csh; echo ./$ndim$SPMHD $rootname$n '>&' $rootname$n.output >> $rootname$n.csh; chmod a+x $rootname$n.csh";
    my $jobid = `cd $rootname; xgrid -hostname cytosine.ex.ac.uk -auth Kerberos -job submit ./$rootname$n.csh`;
##    my $jobid = `cd $rootname; xgrid -hostname cytosine.ex.ac.uk -auth Kerberos -in \$PWD -out \$PWD -se \$PWD/$rootname$n.errors -so \$PWD/$rootname$n.output -job submit ./$ndim$SPMHD $rootname$n`;
    print "$jobid \n";
    ($jobid) = $jobid =~ m/jobIdentifier\s+=\s+(\d+);/; # \s matches spaces (+ = at least one) \d decimals
    print "job id = $jobid \n";
    system "cd $rootname; echo echo getting results of $rootname$n... >> getresults; echo xgrid -hostname cytosine.ex.ac.uk -auth Kerberos -out \$PWD -se \$PWD/$rootname$n.xgriderr -so \$PWD/$rootname$n.xgridout -job results -id $jobid >> getresults$n; echo source getresults$n >> getresults";
    system "cd $rootname; echo echo deleting $rootname$n... >> cleanup; echo xgrid -hostname cytosine.ex.ac.uk -auth Kerberos -job delete -id $jobid >> cleanup";
     system "sleep 2";
}

# call the sametime program to write simultaneous data steps
#system "cd $rootname; ../utils/sametime $rootname $nruns $nsteps";

print "finished $rootname \n";

exit;
