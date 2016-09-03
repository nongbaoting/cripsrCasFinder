#!/usr/bin/perl -w 

use strict;
use FindBin qw/$Bin/;
use File::Spec::Functions qw/rel2abs catfile/;
use File::Basename qw /dirname basename/;
use File::Find qw /find/;
use Cwd qw/cwd/;
use Data::Dumper qw /Dumper/;
use Getopt::Long qw/GetOptions/;

my $cwd = cwd;
my $BinScript = catfile($Bin,"scripts");

my ($genome, $outdir) = @ARGV;
$genome =rel2abs($genome);
die 'argv: my ($genome, $outdir) = @ARGV;' unless @ARGV == 2;

my $result_dir=catfile($cwd ,$outdir);
# `rm -rf $result_dir` if -d $result_dir;
`mkdir -p $result_dir`;

chdir $result_dir;
print "result directory: $result_dir\n";
my $repeat = catfile($result_dir, "repeat.txt"); 
my $out_crispr_around = catfile($result_dir, "crispr_around.fa");

# system "perl $BinScript/1_find_repeat.pl $genome $repeat 15";
# system "perl $BinScript/2_extract_cripsr_around.pl $genome $out_crispr_around $repeat";
`gmsn.pl  --fnn --faa --prok  --output result $out_crispr_around`;
`/data1/exec/ncbi-blast-2.2.28/bin/rpsblast  -query result.faa  -db /data4/melody/source/cdd/cdd/crispr -seg no -evalue 0.01 -outfmt 5 -out cdd_annot.xml`;
`rpsbproc -c /data4/melody/source/cdd/data/rpsbproc.ini -i cdd_annot.xml -o cdd_annot.txt`;
# system "perl $BinScript/4_rps_blast.pl ";
# system "perl $BinScript/5_format_result.pl ";