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

#my ($rh_directories, $dl_jobid);
if (-e "$projectdir/refgenotype/mergedSVs/$mergedirname") {
    print STDERR "Directory $projectdir exists. Please specify a non-existent directory.\n";
    exit;
}
    
# create directories:
    
my $rh_directories = make_directories($projectdir, $mergedirname);
    
my $merge_jobid = launch_targetmerge_swarm($rh_directories, $baitregions);
print "Launched target merge swarm with job id $merge_jobid\n";

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

sub launch_checkdownload {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $dl_jobid = shift;

    my $sample_dir = $rh_dirs->{sampledir};
    my $sample = $sample_dir;
    $sample =~ s/.*\///;

    cp $scriptdir."/sh.checkdownload", $rh_dirs->{"scripts_dir"}."/sh.checkdownload"
        or die "Couldn\'t copy $scriptdir/sh.checkdownload to $rh_dirs->{scripts_dir}!\n";
    my $dependency_string = ($dl_jobid) ? "--dependency=afterok:$dl_jobid" : '';
    system("sbatch $dependency_string --time=48:00:00 --job-name=checkdownload -o $rh_dirs->{log_dir}/\%x_\%j.out -e $rh_dirs->{log_dir}/\%x_\%j.err $rh_dirs->{scripts_dir}/sh.checkdownload $sample $sample_dir > $rh_dirs->{log_dir}/checkdownload.sbatchsubmit.out");

    open CHECK_JOBID, "$rh_dirs->{log_dir}/checkdownload.sbatchsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/checkdownload.sbatchsubmit.out\n";
    my $check_jobid = <CHECK_JOBID>;
    chomp $check_jobid;
    close CHECK_JOBID;

    if ($check_jobid =~ /^\d+$/) {
        return $check_jobid;
    }
    else {
        die "Unable to retrieve job id for check download job in $rh_dirs->{log_dir}/checkdownload.sbatchsubmit.out\n";
    }
}

sub launch_mergeindexminimap {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $dl_jobid = shift;

    my $sample_dir = $rh_dirs->{sampledir};
    my $sample = $sample_dir;
    $sample =~ s/.*\///;

    cp $scriptdir."/sh.mergeandindexminimap", $rh_dirs->{"scripts_dir"}."/sh.mergeandindexminimap";
    my $dependency_string = ($dl_jobid) ? "--dependency=afterok:$dl_jobid" : '';
    system("sbatch $dependency_string --time=24:00:00 $rh_dirs->{scripts_dir}/sh.mergeandindexminimap $sample $sample_dir/minimap2 > $rh_dirs->{log_dir}/read_mm2.sbatchsubmit.out");

    open MERGEINDEX_JOBID, "$rh_dirs->{log_dir}/read_mm2.sbatchsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/read_mm2.sbatchsubmit.out\n";
    my $mergeindex_jobid = <MERGEINDEX_JOBID>;
    chomp $mergeindex_jobid;
    close MERGEINDEX_JOBID;

    if ($mergeindex_jobid =~ /^\d+$/) {
        return $mergeindex_jobid;
    }
    else {
        die "Unable to retrieve job id for merge and index job in $rh_dirs->{log_dir}/read_mm2.sbatchsubmit.out\n";
    }
}

sub launch_createreadtable {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $merge_jobid = shift;
    my $check_jobid = shift;

    my $sample_dir = $rh_dirs->{sampledir};
    my $sample = $sample_dir;
    $sample =~ s/.*\///;

    cp $scriptdir."/sh.create_read_table", $rh_dirs->{"scripts_dir"}."/sh.create_read_table";
    my @dependencies = ();
    foreach my $jobid ($merge_jobid, $check_jobid) {
        if ($jobid) {
            push @dependencies, "afterok:$jobid";
        }
    }
    my $afterok_string = join ",", @dependencies;
    my $dependency_string = (@dependencies) ? "--dependency=$afterok_string" : '';
    system("sbatch $dependency_string --time=24:00:00 --mem=18g $rh_dirs->{scripts_dir}/sh.create_read_table $sample $sample_dir $baitregions > $rh_dirs->{log_dir}/create_reads.sbatchsubmit.out");

    open CREATEREADS_JOBID, "$rh_dirs->{log_dir}/create_reads.sbatchsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/create_reads.sbatchsubmit.out\n";
    my $create_jobid = <CREATEREADS_JOBID>;
    chomp $create_jobid;
    close CREATEREADS_JOBID;

    if ($create_jobid =~ /^\d+$/) {
        return $create_jobid;
    }
    else {
        die "Unable to retrieve job id for create read table job in $rh_dirs->{log_dir}/create_reads.sbatchsubmit.out\n";
    }
}

sub launch_canucorrect {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $create_jobid = shift;

    my $sample_dir = $rh_dirs->{sampledir};
    my $sample = $sample_dir;
    $sample =~ s/.*\///;

    cp $scriptdir."/sh.run_canu_correct", $rh_dirs->{"scripts_dir"}."/sh.run_canu_correct";

    my $swarm_cmds = $rh_dirs->{scripts_dir} . "/sh.launch_canucorrect";

    open HERVS, $baitregions
        or die "Couldn\'t open $baitregions for reading: $!\n";

    my @hervk_names = ();
    while (<HERVS>) {
        chomp;
        my @fields = split /\t/, $_;
        push @hervk_names, $fields[$#fields];
    }
    close HERVS;

    open SWARMCMD, ">$swarm_cmds"
        or die "Couldn\'t open $swarm_cmds for writing: $!\n";

    foreach my $hervk_name (@hervk_names) {
        print SWARMCMD "$rh_dirs->{scripts_dir}/sh.run_canu_correct $sample $hervk_name\n";
    }
    close SWARMCMD;

    my $dependency_string = ($create_jobid) ? "--dependency=afterok:$create_jobid" : '';
    system("swarm $dependency_string -f $swarm_cmds -g 16 -t 4 > $rh_dirs->{log_dir}/canu.swarmsubmit.out");
    
    open CANU_JOBID, "$rh_dirs->{log_dir}/canu.swarmsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/canu.swarmsubmit.out\n";
    my $canu_jobid = <CANU_JOBID>;
    chomp $canu_jobid;
    close CANU_JOBID;

    if ($canu_jobid =~ /^\d+$/) {
        return $canu_jobid;
    }
    else {
        die "Unable to retrieve job id for canu correct job in $rh_dirs->{log_dir}/canu.swarmsubmit.out\n";
    }

}

sub launch_svrefine {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $create_jobid = shift;
    my $canu_jobid = shift;

    my $sample_dir = $rh_dirs->{sampledir};
    my $sample = $sample_dir;
    $sample =~ s/.*\///;

    cp $scriptdir."/sh.aligncorrectedreads", $rh_dirs->{"scripts_dir"}."/sh.aligncorrectedreads";
    
    my $swarm_cmds = $rh_dirs->{scripts_dir} . "/sh.launch_svrefine";

    open HERVS, $baitregions
        or die "Couldn\'t open $baitregions for reading: $!\n";

    my @hervk_names = ();
    while (<HERVS>) {
        chomp;
        my @fields = split /\t/, $_;
        push @hervk_names, $fields[$#fields];
    }
    close HERVS;

    open SWARMCMD, ">$swarm_cmds"
        or die "Couldn\'t open $swarm_cmds for writing: $!\n";

    foreach my $hervk_name (@hervk_names) {
        print SWARMCMD "$rh_dirs->{scripts_dir}/sh.aligncorrectedreads $sample $hervk_name\n";
    }
    close SWARMCMD;

    my @dependencies = ();
    if ($canu_jobid) {
        push @dependencies, "afterany:$canu_jobid";
    }
    if ($create_jobid) {
        push @dependencies, "afterok:$create_jobid";
    }
    my $dependency_string = (@dependencies) ? join ',', @dependencies : "";
    $dependency_string = "--dependency=$dependency_string" if ($dependency_string);

    system("swarm $dependency_string -b 50 -f $rh_dirs->{scripts_dir}/sh.launch_svrefine > $rh_dirs->{log_dir}/launch_svrefine.swarmsubmit.out");
    
    open SVREFINE_JOBID, "$rh_dirs->{log_dir}/launch_svrefine.swarmsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/launch_svrefine.swarmsubmit.out\n";
    my $svrefine_jobid = <SVREFINE_JOBID>;
    chomp $svrefine_jobid;
    close SVREFINE_JOBID;

    if ($svrefine_jobid =~ /^\d+$/) {
        return $svrefine_jobid;
    }
    else {
        die "Unable to retrieve job id for SVrefine job in $rh_dirs->{log_dir}/launch_svrefine.swarmsubmit.out\n";
    }
}

sub launch_svrefinemerge {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $svrefine_jobid = shift;

    my $sample_dir = $rh_dirs->{sampledir};
    my $sample = $sample_dir;
    $sample =~ s/.*\///;

    # this job will merge all of the bam files of corrected reads into one:
    cp $scriptdir."/sh.mergemummerbams", $rh_dirs->{"scripts_dir"}."/sh.mergemummerbams";

    my $dependency_string = ($svrefine_jobid) ? "--dependency=afterok:$svrefine_jobid" : '';
    system("sbatch $dependency_string $rh_dirs->{scripts_dir}/sh.mergemummerbams $sample $sample_dir > $rh_dirs->{log_dir}/mergemummerbams.sbatchsubmit.out");
    
    open SVREFINEMERGE_JOBID, "$rh_dirs->{log_dir}/mergemummerbams.sbatchsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/mergemummerbams.sbatchsubmit.out\n";
    my $svrefinemerge_jobid = <SVREFINEMERGE_JOBID>;
    chomp $svrefinemerge_jobid;
    close SVREFINEMERGE_JOBID;

    if ($svrefinemerge_jobid =~ /^\d+$/) {
        return $svrefinemerge_jobid;
    }
    else {
        die "Unable to retrieve job id for SVrefine job in $rh_dirs->{log_dir}/mergemummerbams.sbatchsubmit.out\n";
    }
}

sub launch_convert {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $svrefinemerge_jobid = shift;

    my $sample_dir = $rh_dirs->{sampledir};
    my $sample = $sample_dir;
    $sample =~ s/.*\///;

    cp $scriptdir."/sh.run_convert_bam", $rh_dirs->{"scripts_dir"}."/sh.run_convert_bam";

    my $dependency_string = ($svrefinemerge_jobid) ? "--dependency=afterok:$svrefinemerge_jobid" : '';
    system("sbatch $dependency_string $rh_dirs->{scripts_dir}/sh.run_convert_bam $sample > $rh_dirs->{log_dir}/convert.sbatchsubmit.out");
    
    open CONVERT_JOBID, "$rh_dirs->{log_dir}/convert.sbatchsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/convert.sbatchsubmit.out\n";
    my $convert_jobid = <CONVERT_JOBID>;
    chomp $convert_jobid;
    close CONVERT_JOBID;

    if ($convert_jobid =~ /^\d+$/) {
        return $convert_jobid;
    }
    else {
        die "Unable to retrieve job id for conversion job in $rh_dirs->{log_dir}/convert.sbatchsubmit.out\n";
    }
}

