#!/usr/local/bin/perl -w
#

use File::Path qw( make_path );
use File::Copy "cp";
use Getopt::Long;

use strict;

our $SCRIPTSDIR='/data/nhansen/HERV_K_catalog/discovery2/scripts';
our $MERGEDDIR='/data/nhansen/HERV_K_catalog/discovery2/refgenotype/mergedSVs';
our $HERVKFILE='/data/nhansen/HERV_K_catalog/discovery2/prep/anchor_baits/all_hervk_baits.bed';

my $Usage = qq!align_alt_allele_reads.pl <--skipmerge> <--skipassemble> <merge directory>\n!;

my ($skipmerge_opt, $skipassemble_opt, $skipalign_opt);
GetOptions( "skipmerge" => \$skipmerge_opt, "skipassemble" => \$skipassemble_opt, "skipalign" => \$skipalign_opt );

$#ARGV==0
    or die "$Usage";

my $mergedir = $ARGV[0];

my $rh_hervk_info = read_hervk_coords($HERVKFILE);

my ($rh_directories, $merge_jobid);
if (!$skipmerge_opt) {
    if (-e "$MERGEDDIR/$mergedir") {
        print STDERR "Directory $MERGEDDIR/$mergedir exists. Please specify a non-existent directory or skip the merge step with --skipmerge option.\n";
        exit;
    }
    
    # create directories:
    
    $rh_directories = make_directories("$MERGEDDIR/$mergedir");
    
    $merge_jobid = launch_merge_swarm($rh_directories, $rh_hervk_info);
    print "Launched merge swarm with job id $merge_jobid\n";
}
else {
    $rh_directories = make_directories("$MERGEDDIR/$mergedir");
}

my ($assembly_jobid, $align_jobid);

if (!$skipassemble_opt) {
    $assembly_jobid = launch_assemble($rh_directories, $rh_hervk_info, $merge_jobid);
    print "Launched assembly job with job id $assembly_jobid\n";
}

$align_jobid = launch_align($rh_directories, $rh_hervk_info, $assembly_jobid);
print "Launched alignment job with job id $align_jobid\n";

########END MAIN#######

sub read_hervk_coords {
    my $hervk_file = shift;

    my %hervk_info = ();
    open HERVKS, $hervk_file
        or die "Couldn\'t open $hervk_file for reading: $!\n";
    while (<HERVKS>) {
        chomp;
        my @fields = split /\t/, $_;
        $hervk_info{$fields[$#fields]} = $_;
    }
    close HERVKS;

    return {%hervk_info};
}

sub make_directories {
    my $maindir = shift;

    my %directories = ();

    ($directories{'maindir'}) = (-e $maindir) ? $maindir :  make_path $maindir;

    ($directories{'scripts_dir'}) = (-e "$maindir/scripts") ? "$maindir/scripts" : make_path "$maindir/scripts";
    ($directories{'log_dir'}) = (-e "$maindir/logs") ? "$maindir/logs" : make_path "$maindir/logs";
    ($directories{'merge_dir'}) = (-e "$maindir/merge") ? "$maindir/merge" : make_path "$maindir/merge";
    ($directories{'assembly_dir'}) = (-e "$maindir/assemble") ? "$maindir/assemble" : make_path "$maindir/assemble";
    ($directories{'assembly_logdir'}) = (-e "$maindir/assemble/logs") ? "$maindir/assemble/logs" : make_path "$maindir/assemble/logs";
    ($directories{'align_dir'}) = (-e "$maindir/align") ? "$maindir/align" : make_path "$maindir/align";

    return {%directories};
}

sub launch_merge_swarm {
    my $rh_dirs = shift;
    my $rh_hervk_info = shift;

    print "MERGEDIR: $rh_dirs->{merge_dir}\n";

    my $merge_dir = $rh_dirs->{merge_dir};

    cp $SCRIPTSDIR."/sh.merge_hervk_svs", $rh_dirs->{"scripts_dir"}."/sh.merge_hervk_svs";
    my $commandfile = $rh_dirs->{"scripts_dir"}."/sh.mergeswarm";
    open COMMANDS, ">$commandfile"
        or die "Couldn\'t open $commandfile for writing: $!\n";

    foreach my $hervkid (keys %{$rh_hervk_info}) {
        print COMMANDS $rh_dirs->{"scripts_dir"}."/sh.merge_hervk_svs $hervkid $merge_dir\n";
    }

    close COMMANDS;

    system("swarm --job-name mergesvs --logdir $rh_dirs->{log_dir} --maxrunning 50 -f $commandfile > $rh_dirs->{log_dir}/mergesvs.swarmsubmit.out");

    open MERGE_JOBID, "$rh_dirs->{log_dir}/mergesvs.swarmsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/mergesvs.swarmsubmit.out\n";
    my $mergejobid = <MERGE_JOBID>;
    chomp $mergejobid;
    close MERGE_JOBID;

    if ($mergejobid =~ /^\d+$/) {
        return $mergejobid;
    }
    else {
        die "Unable to retrieve job id for merge SVs swarm in $rh_dirs->{log_dir}/mergesvs.swarmsubmit.out\n";
    }
}

sub launch_assemble {
    my $rh_dirs = shift;
    my $rh_hervk_info = shift;
    my $merge_jobid = shift;

    print "ASSEMBLYDIR: $rh_dirs->{assembly_dir}\n";

    my $assembly_dir = $rh_dirs->{assembly_dir};

    cp $SCRIPTSDIR."/sh.run_sv_assemble", $rh_dirs->{"scripts_dir"}."/sh.run_sv_assemble";
    my $commandfile = $rh_dirs->{"scripts_dir"}."/sh.assemblyswarm";
    open COMMANDS, ">$commandfile"
        or die "Couldn\'t open $commandfile for writing: $!\n";

    foreach my $hervkid (keys %{$rh_hervk_info}) {
        print COMMANDS $rh_dirs->{"scripts_dir"}."/sh.run_sv_assemble $hervkid $MERGEDDIR/$mergedir\n";
    }

    close COMMANDS;

    my $dependency_string = ($merge_jobid) ? "--dependency=afterany:$merge_jobid" : '';
    system("swarm $dependency_string --job-name svassemble --logdir $rh_dirs->{log_dir} --maxrunning 50 -f $commandfile > $rh_dirs->{log_dir}/svassemble.swarmsubmit.out");

    open ASSEMBLY_JOBID, "$rh_dirs->{log_dir}/svassemble.swarmsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/svassemble.swarmsubmit.out\n";
    my $svassemblejobid = <ASSEMBLY_JOBID>;
    chomp $svassemblejobid;
    close ASSEMBLY_JOBID;

    if ($svassemblejobid =~ /^\d+$/) {
        return $svassemblejobid;
    }
    else {
        die "Unable to retrieve job id for assemble SVs swarm in $rh_dirs->{log_dir}/svassemble.swarmsubmit.out\n";
    }
}

sub launch_align {
    my $rh_dirs = shift;
    my $rh_hervk_info = shift;
    my $assembly_jobid = shift;

    print "ALIGNDIR: $rh_dirs->{align_dir}\n";

    my $align_dir = $rh_dirs->{align_dir};
    my $assembly_dir = $rh_dirs->{assembly_dir};

    opendir ALTS, $assembly_dir
        or die "Couldn\'t open $assembly_dir for reading: $!\n";
    my @alts = grep /ALT/, readdir ALTS;
    closedir ALTS;

    cp $SCRIPTSDIR."/sh.alignaltreads", $rh_dirs->{"scripts_dir"}."/sh.alignaltreads";
    my $commandfile = $rh_dirs->{"scripts_dir"}."/sh.alignswarm";
    open COMMANDS, ">$commandfile"
        or die "Couldn\'t open $commandfile for writing: $!\n";

    foreach my $altdir (@alts) {
        print COMMANDS $rh_dirs->{"scripts_dir"}."/sh.alignaltreads $altdir $MERGEDDIR/$mergedir\n";
    }

    close COMMANDS;

    my $dependency_string = ($assembly_jobid) ? "--dependency=afterany:$assembly_jobid" : '';
    system("swarm $dependency_string --job-name alignaltreads --logdir $rh_dirs->{log_dir} --maxrunning 50 -f $commandfile > $rh_dirs->{log_dir}/alignaltreads.swarmsubmit.out");

    open ALIGN_JOBID, "$rh_dirs->{log_dir}/alignaltreads.swarmsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/alignaltreads.swarmsubmit.out\n";
    my $alignjobid = <ALIGN_JOBID>;
    chomp $alignjobid;
    close ALIGN_JOBID;

    if ($alignjobid =~ /^\d+$/) {
        return $alignjobid;
    }
    else {
        die "Unable to retrieve job id for assemble SVs swarm in $rh_dirs->{log_dir}/alignaltreads.swarmsubmit.out\n";
    }
}
