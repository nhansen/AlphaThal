#!/usr/local/bin/perl -w
#

use File::Path qw( make_path );
use File::Copy "cp";
use Getopt::Long;

use strict;

our $SRAFILE = $ENV{'SRASAMPLEFILE'} || $ENV{'LONGREADTOPDIR'}."/SRA_samples/SRAHumanPacBioWGSRuns.linux.051119.txt';

my $Usage = qq!run_hervk_pipeline.pl <--flankseqs flank_fasta> <--baitregions baits_bed> <--repeatseq repeat_consensus_fasta> <--scripts script_dir> <--skipdownload> <--skipcheck> <--skipmerge> <--skipreadtable> <--skipcanucorrect> <--skipsvrefine> <biosample> <sample directory>\n!;

my $flankseqs = '/data/nhansen/HERV_K_catalog/discovery2/prep/anchor_baits/ref_and_nonref_hervk.flank.fasta';
my $baitregions = '/data/nhansen/HERV_K_catalog/discovery2/prep/anchor_baits/all_hervk_baits.bed';
my $scriptdir='/data/nhansen/HERV_K_catalog/discovery2/scripts';
my $repeatseq = ''; # optional repeat consensus to gather extra reads from

my ($skipdownload_opt, $skipcheck_opt, $skipmerge_opt, $skipreadtable_opt, $skipcanucorrect_opt, $skipsvrefine_opt);
GetOptions("flankseqs=s" => \$flankseqs, "baitregions=s" => \$baitregions, "scripts=s" => \$scriptdir, "repeatseq=s" => \$repeatseq, "skipdownload" => \$skipdownload_opt, "skipcheck" => \$skipcheck_opt, "skipmerge" => \$skipmerge_opt, "skipreadtable" => \$skipreadtable_opt, "skipcanucorrect" => \$skipcanucorrect_opt, "skipsvrefine" => \$skipsvrefine_opt);

$#ARGV==1
    or die "$Usage";

my $biosample = $ARGV[0];
my $sampledir = $ARGV[1];

my ($rh_directories, $dl_jobid);
if (!$skipdownload_opt) {
    if (-e $sampledir) {
        print STDERR "Directory $sampledir exists. Please specify a non-existent directory.\n";
        exit;
    }
    
    # create directories:
    
    $rh_directories = make_directories($sampledir);
    
    $dl_jobid = launch_download_swarm($biosample, $rh_directories);
    print "Launched download swarm with job id $dl_jobid\n";
}
else {
    $rh_directories = make_directories($sampledir);
}

my ($check_jobid, $merge_jobid, $create_jobid, $canu_jobid, $svrefine_jobid, $svrefinemerge_jobid, $convert_jobid);
if (!$skipcheck_opt) {
    $check_jobid = launch_checkdownload($biosample, $rh_directories, $dl_jobid);
    print "Launched check download job with job id $check_jobid\n";
}

if (!$skipmerge_opt) {
    $merge_jobid = launch_mergeindexminimap($biosample, $rh_directories, $dl_jobid);
    print "Launched merge minimap2 job with job id $merge_jobid\n";
}

if (!$skipreadtable_opt) {
    $create_jobid = launch_createreadtable($biosample, $rh_directories, $merge_jobid, $check_jobid);
    print "Launched create read table job with job id $create_jobid\n";
}

if (!$skipcanucorrect_opt) {
    $canu_jobid = launch_canucorrect($biosample, $rh_directories, $create_jobid);
    print "Launched canu correct job with job id $canu_jobid\n";
}

if (!$skipsvrefine_opt) {
    $svrefine_jobid = launch_svrefine($biosample, $rh_directories, $create_jobid, $canu_jobid);
    print "Launched svrefine job with job id $svrefine_jobid\n";
    $svrefinemerge_jobid = launch_svrefinemerge($biosample, $rh_directories, $svrefine_jobid);
    print "Launched mummer BAM merge job with job id $svrefinemerge_jobid\n";
    $convert_jobid = launch_convert($biosample, $rh_directories, $svrefinemerge_jobid);
    print "Launched conversion job with job id $convert_jobid\n";
}

########END MAIN#######

sub make_directories {
    my %directories = ();
    ($directories{'sampledir'}) = (-e $sampledir) ? $sampledir :  make_path $sampledir;

    ($directories{'scripts_dir'}) = (-e "$sampledir/scripts") ? "$sampledir/scripts" : make_path "$sampledir/scripts";
    ($directories{'log_dir'}) = (-e "$sampledir/logs") ? "$sampledir/logs" : make_path "$sampledir/logs";
    ($directories{'data_dir'}) = (-e "$sampledir/SRA_data") ? "$sampledir/SRA_data" : make_path "$sampledir/SRA_data";
    ($directories{'minimap_dir'}) = (-e "$sampledir/minimap2") ? "$sampledir/minimap2" : make_path "$sampledir/minimap2";
    ($directories{'mashmap_dir'}) = (-e "$sampledir/MASHmap") ? "$sampledir/MASHmap" : make_path "$sampledir/MASHmap";
    ($directories{'canu_dir'}) = (-e "$sampledir/canu_correct") ? "$sampledir/canu_correct" : make_path "$sampledir/canu_correct";
    ($directories{'svrefine_dir'}) = (-e "$sampledir/allele_aligns") ? "$sampledir/allele_aligns" : make_path "$sampledir/allele_aligns";
    ($directories{'consensus_dir'}) = (-e "$sampledir/consensus_seqs") ? "$sampledir/consensus_seqs" : make_path "$sampledir/consensus_seqs";
    ($directories{'reports_dir'}) = (-e "$sampledir/Reports") ? "$sampledir/Reports" : make_path "$sampledir/Reports";

    return {%directories};
}

sub launch_download_swarm {
    my $biosample = shift;
    my $rh_dirs = shift;

    print "SAMPLEDIR: $rh_dirs->{sampledir}\n";

    my $sample_dir = $rh_dirs->{sampledir};
    my $accfile = "$sample_dir/SRR_Acc_List.txt";
    my $rh_sample_accessions = retrieve_sample_accessions();

    cp $scriptdir."/sh.download_and_match_to_anchors", $rh_dirs->{"scripts_dir"}."/sh.download_and_match_to_anchors";
    my $commandfile = $rh_dirs->{"scripts_dir"}."/sh.downloadswarm";
    open COMMANDS, ">$commandfile"
        or die "Couldn\'t open $commandfile for writing: $!\n";

    open SRRACC, ">$accfile"
        or die "Couldn\'t open $accfile for writing: $!\n";

    foreach my $srr_acc (keys %{$rh_sample_accessions}) {
        my $platform = $rh_sample_accessions->{$srr_acc}->{'Platform'};
        my $this_biosample = $rh_sample_accessions->{$srr_acc}->{'BioSample'};

        next if ((!$this_biosample) || ($this_biosample ne $biosample));

        print SRRACC "$srr_acc\t$platform\n";
        my $opt_repeat = ($repeatseq) ? " $repeatseq" : "";
        print COMMANDS $rh_dirs->{"scripts_dir"}."/sh.download_and_match_to_anchors $srr_acc $flankseqs $sample_dir \"$platform\"$opt_repeat\n";
    }

    close COMMANDS;
    close SRRACC;

    system("swarm --job-name downloadandmatch --time 32:00:00 --logdir $rh_dirs->{log_dir} --maxrunning 50 -f $commandfile -g 16 > $rh_dirs->{log_dir}/downloadandmatch.swarmsubmit.out");

    open DLAM_JOBID, "$rh_dirs->{log_dir}/downloadandmatch.swarmsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/downloadandmatch.swarmsubmit.out\n";
    my $dl_jobid = <DLAM_JOBID>;
    chomp $dl_jobid;
    close DLAM_JOBID;

    if ($dl_jobid =~ /^\d+$/) {
        return $dl_jobid;
    }
    else {
        die "Unable to retrieve job id for download swarm in $rh_dirs->{log_dir}/downloadandmatch.swarmsubmit.out\n";
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

sub retrieve_sample_accessions {

    open SRA, $SRAFILE
        or die "Couldn\'t open $SRAFILE: $!\n";

    my %srrs = ();
    my @header_fields;
    while (<SRA>) {
        chomp;
        if (/^Run/) {
            @header_fields = split /\t/, $_;
        }
        else {
            my @srr_fields = split /\t/, $_;
            my $srracc = $srr_fields[0];
            $srrs{$srracc} = {};
            for (my $i=0; $i<=$#header_fields; $i++) {
                $srrs{$srracc}->{$header_fields[$i]} = $srr_fields[$i];
            }
        }
    }
    close SRA;

    return {%srrs};
}
