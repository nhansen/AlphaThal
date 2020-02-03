#!/usr/local/bin/perl -w
#

use File::Path qw( make_path );
use File::Copy "cp";
use Getopt::Long;

use strict;

our $TOPDIR = $ENV{'LONGREADTOPDIR'};

if (!$TOPDIR) {
    die "You must set the LONGREADTOPDIR environment variable before running this script!\n";
}

our $SRAFILE = $ENV{'SRASAMPLEFILE'} || $ENV{'LONGREADTOPDIR'}.'/SRASamples/SRAHumanPacBioWGSRuns.linux.051119.txt';

my $Usage = qq!run_longread_pipeline.pl <--flankseqs flank_fasta> <--baitregions baits_bed> <--wideregions wideregions_bed> <--repeatseq repeat_consensus_fasta> <--scripts script_dir> <--localfastadir path_to_fasta_dir> <--reffasta path_to_reference_fasta> <--skipdownload> <--skipcheck> <--skipcanucorrect> <--skipsvrefine> <biosample> <sample directory>\n!;

my $flankseqs;
my $baitregions;
my $scriptdir="$TOPDIR/scripts";
my $reffasta;
my $wideregions = "$TOPDIR/prep/ref_and_nonref.wideregions.bed";
my $repeatseq = ''; # optional repeat consensus to gather extra reads from
my $localfastadir_opt; # optional local directory with long read fasta files to screen
my $platform_opt; # option to specify the platform of the sequence

my ($skipdownload_opt, $skipcheck_opt, $skipcanucorrect_opt, $skipsvrefine_opt);
GetOptions("flankseqs=s" => \$flankseqs, "baitregions=s" => \$baitregions, "wideregions=s" => \$wideregions, "scripts=s" => \$scriptdir, "repeatseq=s" => \$repeatseq, "skipdownload" => \$skipdownload_opt, "skipcheck" => \$skipcheck_opt, "skipcanucorrect" => \$skipcanucorrect_opt, "skipsvrefine" => \$skipsvrefine_opt, "localfastadir=s" => \$localfastadir_opt, "reffasta=s" => \$reffasta, "platform=s" => \$platform_opt);

$#ARGV==1
    or die "$Usage";

my $biosample = $ARGV[0];
my $sampledir = $ARGV[1];

if (!$baitregions || !$flankseqs || !$reffasta) {
    die "Must specify a bed file of bait regions with --baitregions, an indexed fasta file of flanking seqs with --flankseqs, and a reference fasta file with --reffasta!\n$Usage";
}

my ($rh_directories, $dl_jobid);

if ($localfastadir_opt) { # skip download and launch just the MASHmap/minimap jobs on specified fasta files

    # create directories:
    $rh_directories = make_directories($sampledir);

    if (!$skipdownload_opt) {

        $dl_jobid = launch_match_swarm($biosample, $rh_directories, $localfastadir_opt, $reffasta, $baitregions, $platform_opt);
        print "Launched match swarm with job id $dl_jobid\n";
    }
}
elsif (!$skipdownload_opt) {
    if (-e $sampledir) {
        print STDERR "Directory $sampledir exists. Please specify a non-existent directory.\n";
        exit;
    }
    
    # create directories:
    
    $rh_directories = make_directories($sampledir);
    
    $dl_jobid = launch_download_swarm($biosample, $rh_directories, $reffasta);
    print "Launched download swarm with job id $dl_jobid\n";
}
else { # probably have alignment results already
    $rh_directories = make_directories($sampledir);
}

my ($check_jobid, $canu_jobid, $svrefine_jobid, $mummerbammerge_jobid, $convert_jobid);
if (!$skipcheck_opt) {
    $check_jobid = launch_checkdownload($biosample, $rh_directories, $dl_jobid);
    print "Launched check download job with job id $check_jobid\n";
}

if (!$skipcanucorrect_opt) {
    $canu_jobid = launch_canucorrect($biosample, $rh_directories, $check_jobid);
    print "Launched canu correct job with job id $canu_jobid\n";
}

if (!$skipsvrefine_opt) {
    $svrefine_jobid = launch_svrefine($biosample, $rh_directories, $canu_jobid);
    print "Launched svrefine job with job id $svrefine_jobid\n";
    $mummerbammerge_jobid = launch_mummerbammerge($biosample, $rh_directories, $svrefine_jobid);
    print "Launched mummer BAM merge job with job id $mummerbammerge_jobid\n";
    $convert_jobid = launch_convert($biosample, $rh_directories, $wideregions, $reffasta, $mummerbammerge_jobid);
    print "Launched conversion job with job id $convert_jobid\n";
}

########END MAIN#######

sub make_directories {
    my %directories = ();
    ($directories{'sampledir'}) = (-e $sampledir) ? $sampledir :  make_path $sampledir;

    ($directories{'scripts_dir'}) = (-e "$sampledir/scripts") ? "$sampledir/scripts" : make_path "$sampledir/scripts";
    ($directories{'log_dir'}) = (-e "$sampledir/logs") ? "$sampledir/logs" : make_path "$sampledir/logs";
    ($directories{'data_dir'}) = (-e "$sampledir/read_data") ? "$sampledir/read_data" : make_path "$sampledir/read_data";
    ($directories{'minimap_dir'}) = (-e "$sampledir/minimap2") ? "$sampledir/minimap2" : make_path "$sampledir/minimap2";
    ($directories{'mashmap_dir'}) = (-e "$sampledir/MASHmap") ? "$sampledir/MASHmap" : make_path "$sampledir/MASHmap";
    ($directories{'canu_dir'}) = (-e "$sampledir/canu_correct") ? "$sampledir/canu_correct" : make_path "$sampledir/canu_correct";
    ($directories{'readtable_dir'}) = (-e "$sampledir/read_table") ? "$sampledir/read_table" : make_path "$sampledir/read_table";
    ($directories{'svrefine_dir'}) = (-e "$sampledir/allele_aligns") ? "$sampledir/allele_aligns" : make_path "$sampledir/allele_aligns";
    #($directories{'consensus_dir'}) = (-e "$sampledir/consensus_seqs") ? "$sampledir/consensus_seqs" : make_path "$sampledir/consensus_seqs";
    #($directories{'reports_dir'}) = (-e "$sampledir/Reports") ? "$sampledir/Reports" : make_path "$sampledir/Reports";

    return {%directories};
}

sub launch_download_swarm {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $ref_fasta = shift;

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
        print COMMANDS $rh_dirs->{"scripts_dir"}."/sh.download_and_match_to_anchors $srr_acc $flankseqs $ref_fasta $sample_dir \"$platform\"$opt_repeat\n";
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

sub launch_match_swarm {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $fastadir = shift;
    my $reffasta = shift;
    my $baits = shift;
    my $platform = shift;

    print "SAMPLEDIR: $rh_dirs->{sampledir}\n";

    my $sample_dir = $rh_dirs->{sampledir};
    my $fileoffiles = "$sample_dir/read_fasta_list.txt"; # create a file with fasta locations

    cp $scriptdir."/sh.match_to_anchors", $rh_dirs->{"scripts_dir"}."/sh.match_to_anchors";
    my $commandfile = $rh_dirs->{"scripts_dir"}."/sh.matchswarm";
    open COMMANDS, ">$commandfile"
        or die "Couldn\'t open $commandfile for writing: $!\n";

    open FASTAFILES, ">$fileoffiles"
        or die "Couldn\'t open $fileoffiles for writing: $!\n";

    opendir FASTAS, $fastadir
        or die "Couldn\'t open directory $fastadir for reading: $!\n";

    my @fastafiles = grep /\.fa(sta){0,1}(\.gz){0,1}$/, readdir FASTAS;

    closedir FASTAS;

    foreach my $fastafile (@fastafiles) {
        my $platform = $platform_opt || 'Unknown';

        print FASTAFILES "$fastadir/$fastafile\t$platform\n";

        symlink "$fastadir/$fastafile", "$rh_dirs->{data_dir}/$fastafile"
            or die "Couldn\'t create symbolic link to $fastadir/$fastafile in $rh_dirs->{data_dir}: $!\n";
        my $opt_repeat = ($repeatseq) ? " $repeatseq" : "";
        print COMMANDS $rh_dirs->{"scripts_dir"}."/sh.match_to_anchors $fastafile $flankseqs $reffasta $baits $sample_dir $opt_repeat\n";
    }

    close COMMANDS;
    close FASTAFILES;

    system("swarm --job-name matchtoanchors --time 32:00:00 --logdir $rh_dirs->{log_dir} -f $commandfile -g 16 > $rh_dirs->{log_dir}/matchtoanchors.swarmsubmit.out");

    open DLAM_JOBID, "$rh_dirs->{log_dir}/matchtoanchors.swarmsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/matchtoanchors.swarmsubmit.out\n";
    my $dl_jobid = <DLAM_JOBID>;
    chomp $dl_jobid;
    close DLAM_JOBID;

    if ($dl_jobid =~ /^\d+$/) {
        return $dl_jobid;
    }
    else {
        die "Unable to retrieve job id for match swarm in $rh_dirs->{log_dir}/matchtoanchors.swarmsubmit.out\n";
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

sub launch_canucorrect {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $check_jobid = shift;

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

    my $dependency_string = ($check_jobid) ? "--dependency=afterok:$check_jobid" : '';
    system("swarm $dependency_string -f $swarm_cmds --logdir $rh_dirs->{log_dir} -b 5 -g 16 -t 4 > $rh_dirs->{log_dir}/canu.swarmsubmit.out");
    
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
    my $dependency_string = (@dependencies) ? join ',', @dependencies : "";
    $dependency_string = "--dependency=$dependency_string" if ($dependency_string);

    system("swarm $dependency_string -b 50 --logdir $rh_dirs->{log_dir} -f $rh_dirs->{scripts_dir}/sh.launch_svrefine > $rh_dirs->{log_dir}/launch_svrefine.swarmsubmit.out");
    
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

sub launch_mummerbammerge {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $svrefine_jobid = shift;

    my $sample_dir = $rh_dirs->{sampledir};
    my $sample = $sample_dir;
    $sample =~ s/.*\///;

    # this job will merge all of the bam files of corrected reads into one:
 
    opendir ALIGNS, $rh_dirs->{svrefine_dir}
        or die "Couldn\'t open $rh_dirs->{svrefine_dir} for reading: $!\n";

    my @mummerbams = grep /\.fullheader\.sort\.bam$/, readdir ALIGNS;
    closedir ALIGNS;

    my $numbams = @mummerbams;
    if ($numbams <= 1000) {
        cp $scriptdir."/sh.mergemummerbams", $rh_dirs->{"scripts_dir"}."/sh.mergemummerbams";
    }
    else { # write our own merge script for a multi-pass merge:
        my $mergescript = "$rh_dirs->{'scripts_dir'}/sh.mergemummerbams";
        open MERGE, ">$mergescript"
            or die "Couldn\'t open $mergescript for writing: $!\n";
        print MERGE "#!/bin/bash\n\nmodule load samtools\n\n";
        print MERGE "export SAMPLE=\$1\nexport SAMPLEDIR=\$2\n\ncd \$SAMPLEDIR/allele_aligns\n";

        my @outputfiles = ();
        my $outputfilenum = 1;
        while (@mummerbams) {
           my $thisoutput = "$sample.genome.sort.$outputfilenum.bam";
           $outputfilenum++;
           print MERGE "samtools merge -f $thisoutput ";
           for (my $i=1; $i<= 1000; $i++) {
               my $inputfile = shift @mummerbams;
               if ($inputfile) {
                   print MERGE " $inputfile";
               }
           }
           print MERGE "\n";
           push @outputfiles, $thisoutput;
        }

        my $output_bam_string = join ' ', @outputfiles;
        print MERGE "\nsamtools merge -f \$SAMPLE.genome.bam $output_bam_string\n";
        print MERGE "samtools sort \$SAMPLE.genome.bam -o \$SAMPLE.genome.sort.bam\n";
        print MERGE "samtools index \$SAMPLE.genome.sort.bam\n";
        close MERGE;
    }
    
    my $dependency_string = ($svrefine_jobid) ? "--dependency=afterok:$svrefine_jobid" : '';
    system("sbatch $dependency_string -o $rh_dirs->{log_dir}/\%x_\%j.out -e $rh_dirs->{log_dir}/\%x_\%j.err --job-name=mergesvrefine $rh_dirs->{scripts_dir}/sh.mergemummerbams $sample $sample_dir > $rh_dirs->{log_dir}/mergemummerbams.sbatchsubmit.out");
    
    open SVREFINEMERGE_JOBID, "$rh_dirs->{log_dir}/mergemummerbams.sbatchsubmit.out"
        or die "Couldn\'t open $rh_dirs->{log_dir}/mergemummerbams.sbatchsubmit.out\n";
    my $mummerbammerge_jobid = <SVREFINEMERGE_JOBID>;
    chomp $mummerbammerge_jobid;
    close SVREFINEMERGE_JOBID;

    if ($mummerbammerge_jobid =~ /^\d+$/) {
        return $mummerbammerge_jobid;
    }
    else {
        die "Unable to retrieve job id for SVrefine job in $rh_dirs->{log_dir}/mergemummerbams.sbatchsubmit.out\n";
    }
}

sub launch_convert {
    my $biosample = shift;
    my $rh_dirs = shift;
    my $widebedfile = shift;
    my $reffasta = shift;
    my $mummerbammerge_jobid = shift;

    my $sample_dir = $rh_dirs->{sampledir};
    my $sample = $sample_dir;
    $sample =~ s/.*\///;

    cp $scriptdir."/sh.run_convert_bam", $rh_dirs->{"scripts_dir"}."/sh.run_convert_bam";

    my $dependency_string = ($mummerbammerge_jobid) ? "--dependency=afterok:$mummerbammerge_jobid" : '';
    system("sbatch $dependency_string --job-name=convertbam -o $rh_dirs->{log_dir}/\%x_\%j.out -e $rh_dirs->{log_dir}/\%x_\%j.err $rh_dirs->{scripts_dir}/sh.run_convert_bam $sample $widebedfile $reffasta > $rh_dirs->{log_dir}/convert.sbatchsubmit.out");
    
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
