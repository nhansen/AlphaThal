#!/usr/local/bin/perl -w
#
use strict;

use FileHandle;

my $Usage = "multi_merge_bam.pl <file_of_bam_files> <# to merge at a time> <output file>\n";

$#ARGV==2
    or die "$Usage";

my $file_of_files = $ARGV[0];
my $no_to_merge = $ARGV[1];
my $final_output_bam = $ARGV[2];

my $fofstring = ($file_of_files =~ /\.gz$/) ? "gunzip -c $file_of_files | " : $file_of_files;
my $fof_fh = FileHandle->new($fofstring);

my @filenames = ();
while (<$fof_fh>) {
    next if (/^#/);

    if (/^\s*(\S+)\s*/) {
        push @filenames, $1;
    }
}
close $fof_fh;

my @final_files_to_merge = ();
my @smallgroup_to_merge = ();

while (@filenames) {
    push @smallgroup_to_merge, shift @filenames;
    my $smallgroup_size = @smallgroup_to_merge;
    if (($smallgroup_size == $no_to_merge) || (!@filenames)) {
        if ($smallgroup_size == 1) { # last file
            push @final_files_to_merge, $smallgroup_to_merge[0];
        }
        else {
            my $inputstring = join ' ', @smallgroup_to_merge;
            my $mergenum = @final_files_to_merge;
            $mergenum++;
            my $outputbam = 'smallmerge.'.$mergenum.'.bam';
            system("samtools merge -o $outputbam $inputstring")==0
                or die "Couldn\'t merge: samtools merge -o $outputbam $inputstring\n";
            push @final_files_to_merge, $outputbam;
        }
        @smallgroup_to_merge = ();
    }
}

# merge the final set:

my $inputstring = join ' ', @final_files_to_merge;
system("samtools merge -o $final_output_bam $inputstring")==0
    or die "Unable to merge merged files $inputstring into final bam $final_output_bam\n";

