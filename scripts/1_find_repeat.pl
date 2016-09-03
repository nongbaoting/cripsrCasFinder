#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;
use Cwd;
use Parallel::ForkManager;

#my $seed=15;
my ($fasta,$outfile,$seed) = @ARGV;
my $in=Bio::SeqIO->new(-file=>"$fasta",
						-format=>"fasta"
						);

open OUT,">",$outfile;

my ( %results);

while(my $seq_ojb=$in->next_seq){
		my $len=$seq_ojb->length;
		my $seq = $seq_ojb->seq;
		my $seq_id = $seq_ojb->id;
		$seq=uc($seq);
		
		# print "$id\t$len\n";
		my %patterns;
		if($len > 50000){
			 # print "\n$seq_id\t$len\n";
			
			my @seqs = divide($seq,30000,10000);
			my $cpu = 10;
			my $pm = Parallel::ForkManager->new($cpu);
			# data structure retrieval and handling
			$pm -> run_on_finish ( # called BEFORE the first call to start()
				sub {
					my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
  
					# retrieve data structure from child
					if (defined($data_structure_reference)) {  # children are not forced to send anything
						my %pattern_child = %{$data_structure_reference};  # child passed a string reference
						%patterns =(%patterns,%pattern_child) ;
					} else {  # problems occuring during storage or retrieval will throw a warning
					print qq|No message received from child process $pid!\n|;
					}
					}
				);
			foreach my $subseq (@seqs){
				
				my $pid = $pm->start and next;
				my %pattern_child=find_pattern($subseq,$seed);
				# foreach my $key (keys %pattern_child){print "$key !!!!!!!!!!\n";}			
				$pm->finish(0,\%pattern_child);
				
			}	
			$pm->wait_all_children; 
		}else{
			 %patterns = find_pattern($seq, $seed);
		}
		
		
		foreach my $pat (keys %patterns){
			my $length_new = length($pat);
			# print "$pat !!!\n";
			my @pos;
			while($seq=~m/$pat/g){
				my $pos=pos($seq);
				push(@pos,$pos-$length_new+1);
			}
			my @dis;
			my $k=@pos-1;
			@pos = sort {$a <=>$b}@pos;
			my $index;
			for(my $i=0;$i<$k;$i++){
				push(@dis,$pos[$i+1]-$pos[$i]-$length_new);
				$index.=substr($seq,($pos[$i]+$length_new),($pos[$i+1]-$pos[$i]-$length_new))."\t";
				}
				my @sorted=sort {$a<=>$b}@dis;
				my $avg=average(@sorted); 
			my $temp=@dis;
			$temp=int(($temp+1)/2);#print $temp,"\n";
			my $media=$sorted[$temp]; 
			if($media<=100 or $avg<=100){
					 my $key=$pat."\t"."$seq_id";
					 
				 if($pat=~/TCTAC/ and $pat=~/GTAGA/){
					$results{$key} ="1\t$seq_id\t". "@pos"."\t"."@dis"."\t$length_new\t". ($#pos+1) ."\t$pat\t$index\n";}	
						else{
					$results{$key} ="0\t$seq_id\t". "@pos"."\t"."@dis"."\t$length_new\t". ($#pos+1) ."\t$pat\t$index\n";}
					
				}else{
					print "$pat:	Distance too big\n";
				}
		
		
		}
				
		# print join "\t",$id, @pos,"\n";		
}

foreach my $keys (keys %results){
	print OUT "$results{$keys}";
}
print "$fasta\tfinished!\n"; 

sub average {
my (@num) = @_;
my $num = scalar @num;
my $total;
foreach (0..$#num) {

$total += $num[$_];
}
return ($total/$num);
}

sub divide {
		my($seq,$segment,$overlap) = @_;
		my $len = length($seq);
		my $start = 0;
		my @seqs;
		while(1){ 
			my $subseq = substr($seq,$start ,($segment + $overlap));
			push @seqs,$subseq;
			last unless ($start + $segment-$overlap) <$len ;
			$start =$start + $segment-$overlap;
		}
		return @seqs;
}


sub find_pattern {
	my ($seq,$seed) = @_;
	
	my $l = length($seq);
	my $startPos = 0;
	# print "\t$l\n";

	my %pattern;
	while(($startPos + $seed) < $l){
			my $seq15 = substr($seq,$startPos,$seed);
			my @counts = $seq=~ /$seq15/g;
		
			if( @counts >= 3){	
	
				
				my $flag=1;
				while($flag){
				 my $pat = substr($seq,$startPos,($seed+$flag));
				 my $length_new = length($pat);
				 my @counts_n = $seq=~ /$pat/g;
				 # print "$pat !!!\n";
				 #$start = $end - $len + 1;
				# print "$#counts $#counts_n";
				 if($#counts <= $#counts_n){
					$flag++;
					$flag = 0 if ($startPos + $seed + $flag  ) >= $l;
				#	print $pat ,"\n";
				 }else{
					$flag--;
					my $pat = substr($seq,$startPos,($seed+$flag));				
					$startPos =$startPos + $seed+$flag;
					# print "$pat\n";
					$flag=0;
					$pattern{$pat} =1;
					}
				 
				 }
			 }
			$startPos += 1;
	}

	return %pattern;
}