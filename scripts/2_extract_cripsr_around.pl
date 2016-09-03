#!/usr/bin/perl -w 

use strict;
use FindBin qw/$Bin/;
use File::Spec::Functions qw/rel2abs catfile/;
use File::Basename qw /dirname basename/;
use File::Find qw /find/;
use Cwd qw/cwd/;
use Bio::SeqIO;
use Data::Dumper qw /Dumper/;
use Getopt::Long qw/GetOptions/;

my $cwd = cwd;
my ($genome ,$out_crispr_around, $repInfo) = @ARGV;

open my $fh_r, "<", $repInfo or die "can not open file $repInfo\n";
open my $fh_o, ">", $out_crispr_around or die "can not open file $out_crispr_around\n";

my (%fastaHash);

my $in = Bio::SeqIO->new(
					-file => $genome,
					-format => "fasta"	
						);
						
while(my $seq_ojb = $in->next_seq){
	my $id = $seq_ojb->id;
	my $seq = $seq_ojb->seq;
	my $length = $seq_ojb->length;
	$fastaHash{$id} = [$length,$seq];
}

while(<$fh_r>){
	chomp;
	my @xl = split(/\t/,$_);
	my $id = $xl[1];
	my @pos = split(/\s/,$xl[2]);
	next unless defined $fastaHash{$id};
	my $repeat = $xl[6];
	my $crispr_start = $pos[0];
	my $crispr_end = $pos[-1];
	my @values = @{$fastaHash{$id}};
	
	my($up_len, $down_len,$sublen);
	
	if($crispr_start<1000){
		$up_len = -1;
	}elsif($crispr_start <20000){
		$up_len = $crispr_start;
	}else{
		$up_len = 20000;
	}
	
	if( ($values[0] - $crispr_end) <1000){
		$down_len = -1;
	}elsif(($values[0] - $crispr_end) < 20000){
		$down_len = $values[0];
	}else{
		$down_len =  20000; 
	}
	
	unless($up_len == -1){
		my $seq_up = substr($values[1],$up_len,$crispr_start);
		$sublen = length($seq_up);
		print $fh_o "\>$id\_up\_$repeat\_$sublen\n$seq_up\n";
	}
	
	unless($down_len == -1){
		my $seq_down = substr($values[1],$crispr_end, $down_len);
		$sublen = length($seq_down); 
		print $fh_o "\>$id\_down\_$repeat\_$sublen\n$seq_down\n";
	}
}

