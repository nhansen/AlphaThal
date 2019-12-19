#!/usr/local/bin/perl -w
#

use File::Path qw( make_path );
use File::Copy "cp";
use Getopt::Long;

use strict;

my $Usage = qq!merge_project_alt_alleles.pl <--baitregions baits_bed> <--scripts script_dir> <project directory> <merge directory name>\n!;

my $baitregions = '/data/nhansen/HERV_K_catalog/SVA_discovery/prep/anchor_baits/all_hervk_baits.bed';
my $scriptdir='/data/nhansen/HERV_K_catalog/SVA_discovery/scripts';
my $repeatseq = ''; # optional repeat consensus to gather extra reads from

GetOptions( "baitregions=s" => \$baitregions, "scripts=s" => \$scriptdir, "repeatseq=s" => \$repeatseq );

$#ARGV==1
    or die "$Usage";

my $projectdir = $ARGV[0];
my $mergedirname = $ARGV[1];

if (-e "$projectdir/refgenotype/mergedSVs/$mergedirname") {
    print STDERR "Directory $projectdir exists. Please specify a non-existent directory.\n";
    exit;
}
    
# create directories:
    
my $rh_directories = make_directories($projectdir, $mergedirname);
    
my $merge_jobid = launch_targetmerge_swarm($rh_directories, $baitregions);
print "Launched target merge swarm with job id $merge_jobid\n";
my $waa_jobid = launch_writealleles_swarm($rh_directories, $mergedirname, $baitregions, $merge_jobid);
print "Launched alt allele fasta swarm with job id $waa_jobid\n";

########END MAIN#######

sub make_directories {
    my $projectdir = shift;
    my $mergedirname = shift;

    my %directories = ();
    $directories{'projectdir'} = $projectdir;

    my $mergedir = "$projectdir/refgenotype/mergedSVs/$mergedirname";

    ($directories{'mergedir'}) = (-e "$mergedir") ? $mergedir :  make_path $mergedir;
    ($directories{'scriptsdir'}) = (-e "$mergedir/scripts") ? "$mergedir/scripts" : make_path "$mergedir/scripts";
    ($directories{'targetvcfdir'}) = (-e "$mergedir/targetvcfs") ? "$mergedir/targetvcfs" : make_path "$mergedir/targetvcfs";
    ($directories{'altalleledir'}) = (-e "$mergedir/altalleleinfo") ? "$mergedir/altalleleinfo" : make_path "$mergedir/altalleleinfo";
    ($directories{'logdir'}) = (-e "$mergedir/logs") ? "$mergedir/logs" : make_path "$mergedir/logs";

    return {%directories};
}

sub launch_targetmerge_swarm {
    my $rh_dirs = shift;
    my $baitregion_file = shift;

    my $mergedir = $rh_dirs->{mergedir};
    print "MERGEDIR: $mergedir\n";

    my $targetvcf_dir = $rh_dirs->{targetvcfdir};

    cp $scriptdir."/sh.merge_target_svs", $rh_dirs->{"scriptsdir"}."/sh.merge_target_svs";
    my $commandfile = $rh_dirs->{"scriptsdir"}."/sh.mergetargetswarm";
    open COMMANDS, ">$commandfile"
        or die "Couldn\'t open $commandfile for writing: $!\n";

    my $projectdir = $rh_dirs->{projectdir};
    open BAITS, "$baitregion_file"
        or die "Couldn\'t open $baitregion_file for reading: $!\n";

    while (<BAITS>) {
        chomp;
        my @fields = split /\t/, $_;
        my $targetname = $fields[$#fields];

        print COMMANDS $rh_dirs->{"scriptsdir"}."/sh.merge_target_svs $targetname $targetvcf_dir\n";
    }

    close COMMANDS;
    close BAITS;

    system("swarm --job-name mergetargetvcfs --time 2:00:00 --logdir $rh_dirs->{logdir} -b 30 -f $commandfile > $rh_dirs->{logdir}/mergetargetvcfs.swarmsubmit.out");

    open MTV_JOBID, "$rh_dirs->{logdir}/mergetargetvcfs.swarmsubmit.out"
        or die "Couldn\'t open $rh_dirs->{logdir}/mergetargetvcfs.swarmsubmit.out\n";
    my $mtv_jobid = <MTV_JOBID>;
    chomp $mtv_jobid;
    close MTV_JOBID;

    if ($mtv_jobid =~ /^\d+$/) {
        return $mtv_jobid;
    }
    else {
        die "Unable to retrieve job id for merge target vcf swarm in $rh_dirs->{logdir}/mergetargetvcfs.swarmsubmit.out\n";
    }
}

# Use merged VCFs to write fasta files of alt alleles to align reads to:
#
sub launch_writealleles_swarm {
    my $rh_dirs = shift;
    my $mergedirname = shift;
    my $baitregion_file = shift;
    my $mtv_jobid = shift;

    my $mergedir = $rh_dirs->{mergedir};
    print "MERGEDIR: $mergedir\n";

    cp $scriptdir."/sh.write_alt_alleles", $rh_dirs->{"scriptsdir"}."/sh.write_alt_alleles";
    my $commandfile = $rh_dirs->{"scriptsdir"}."/sh.altalleleswarm";
    open COMMANDS, ">$commandfile"
        or die "Couldn\'t open $commandfile for writing: $!\n";

    my $projectdir = $rh_dirs->{projectdir};
    open BAITS, "$baitregion_file"
        or die "Couldn\'t open $baitregion_file for reading: $!\n";

    while (<BAITS>) {
        chomp;
        my @fields = split /\t/, $_;
        my $targetname = $fields[$#fields];

        print COMMANDS $rh_dirs->{"scriptsdir"}."/sh.write_alt_alleles $targetname $mergedirname\n";
    }

    close COMMANDS;
    close BAITS;

    system("swarm --dependency=afterok:$mtv_jobid --job-name writealtalleles --time 2:00:00 --logdir $rh_dirs->{logdir} -b 50 -f $commandfile > $rh_dirs->{logdir}/writealtalleles.swarmsubmit.out");

    open WAA_JOBID, "$rh_dirs->{logdir}/writealtalleles.swarmsubmit.out"
        or die "Couldn\'t open $rh_dirs->{logdir}/writealtalleles.swarmsubmit.out\n";
    my $waa_jobid = <WAA_JOBID>;
    chomp $waa_jobid;
    close WAA_JOBID;

    if ($waa_jobid =~ /^\d+$/) {
        return $waa_jobid;
    }
    else {
        die "Unable to retrieve job id for write alt allele swarm in $rh_dirs->{logdir}/writealtalleles.swarmsubmit.out\n";
    }
}
