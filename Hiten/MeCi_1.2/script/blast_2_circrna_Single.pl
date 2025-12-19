#!/usr/bin/env perl
use strict;
use warnings;

#Author: guantao.zheng
#Create: 20210210

my $indir   = $ARGV[0];
my $outdir  = $ARGV[1];
my $prefix  = $ARGV[2];
my $overlap = $ARGV[3];
my $gap     = $ARGV[4];
my $fasta   = $ARGV[5];
print scalar(@ARGV),"\n";
if(scalar(@ARGV) < 6){
	die "usage: $0 indir outdir prefix overlap gap fasta !";
}

my $max_len = 10000;


my %info_reads;
my %info_circs;
open FOUT,">$outdir/poly(A)_junction_reads.bed";
foreach my $file(glob("$indir/*.txt")){
	open  FIN,$file;
	while(my $line = <FIN>){
		chomp($line);
		my ($read,$strand,$q1_s,$q1_e,$t1_s,$t1_e,$q2_s,$q2_e,$t2_s,$t2_e) = split /\t/,$line;
		my $diff = $q1_e - $q2_s + 1;
		next if $diff < -$gap;
		next if $diff > $overlap;
		if($strand eq "+"){
			next if ($t1_e - $t2_s) > $max_len;
			next if $t1_s <= $t2_e;
			# next if ($t1_s <= $t2_s || $t1_e <= $t2_e);
			# next if $t1_s <= ($t2_e - 3);

			my $start_a = $t2_s;
			my $start_b = $t2_s + $diff;
			my $end_a   = $t1_e;
			my $end_b   = $t1_e;
			my $length  = $end_a - $start_a + 1;

			my $id = "$prefix|$start_a|$end_a|$diff|+";
			$info_circs{$id}{start}  = $start_a;
			$info_circs{$id}{end}    = $end_a;
			$info_circs{$id}{strand} = $strand;
			$info_circs{$id}{reads}{$read}++;
			$info_reads{$read}{circrna}{$id} = $line;
			if($diff < 0){
				print FOUT $read, "\t", $q1_e, "\t", $q2_s - 1, "\t", $id, "\t", $diff, "\t", "+", "\n";
			}
			
		}else{
			next if ($t2_e - $t1_s) > $max_len;
			next if $t2_s <= $t1_e;
			# next if ($t2_s <= $t1_s || $t2_e <= $t1_e);
			# next if $t2_s <= ($t1_e - 3);

			my $start_a = $t2_e;
			my $start_b = $t1_s + $diff;
			my $end_a   = $t1_s;
			my $end_b   = $t2_e;
			my $length = $end_a - $start_a + 1;

			my $id = "$prefix|$start_a|$end_a|$diff|-";
			$info_circs{$id}{start}  = $start_a;
			$info_circs{$id}{end}    = $end_a;
			$info_circs{$id}{strand} = $strand;
			$info_circs{$id}{reads}{$read}++;
			$info_reads{$read}{circrna}{$id} = $line;
			if($diff < 0){
				print FOUT $read, "\t", $q1_e, "\t", $q2_s - 1, "\t", $id, "\t", $diff, "\t", "+", "\n";
			}
		}
	}
	close FIN;
}
close FOUT;

system("bedtools getfasta -fi '$fasta' -fo '$outdir/poly(A)_junction_reads.fa' -bed '$outdir/poly(A)_junction_reads.bed' -name");

my %info_insert;
open  FIN,"$outdir/poly(A)_junction_reads.fa";
my $read_ID = "";
my $circ_ID = "";
while(<FIN>){
	chomp;
	if(/^>(.+)/){
		$read_ID = "";
		$circ_ID = "";
		my $ttt = $1;
		my @p = split /:/,$ttt;
		my $circ_id = shift @p;
		shift @p;
		my $tmp = pop @p;
		$read_ID = join(":",@p);
		my ($s,$e) = split /-/,$tmp;
		# $insert_ID = $circ_id."__".$s."__".$e;
		$circ_ID = $circ_id;
		$info_insert{$read_ID}{$circ_ID}{gene} = $circ_id;
	}else{
		my $seq = $_;
		my $len = length($seq);
		my $num_A = $seq =~ tr/A//;
		my $num_T = $seq =~ tr/T//;
		my $num_G = $seq =~ tr/G//;
		my $num_C = $seq =~ tr/C//;
		$info_insert{$read_ID}{$circ_ID}{sequence} = $seq;
		$info_insert{$read_ID}{$circ_ID}{length}   = $len;
		$info_insert{$read_ID}{$circ_ID}{percent_A}  = sprintf("%.2f",$num_A / $len * 100);
		$info_insert{$read_ID}{$circ_ID}{percent_T}  = sprintf("%.2f",$num_T / $len * 100);
		$info_insert{$read_ID}{$circ_ID}{percent_G}  = sprintf("%.2f",$num_G / $len * 100);
		$info_insert{$read_ID}{$circ_ID}{percent_C}  = sprintf("%.2f",$num_C / $len * 100);
		# print $insert_ID,"\n";
		# print $read_ID,"\t",$insert_ID,"\t",$ratio,"\n";
	}

}
close FIN;

foreach my $read(sort keys %info_reads){
	if(scalar keys %{$info_reads{$read}{circrna}} > 1){
		$info_reads{$read}{type} = "multiple";
	}else{
		$info_reads{$read}{type} = "unique";
	}
}

open FOUT,">$outdir/circrna_details.xls";
print FOUT "circrna_id\tchr\tstart\tend\tstrand\tcount\treads\n";
foreach my $id(sort keys %info_circs){
	print FOUT $id,"\t",$prefix,"\t",$info_circs{$id}{start},"\t",$info_circs{$id}{end},"\t",$info_circs{$id}{strand},"\t",scalar(keys %{$info_circs{$id}{reads}}),"\t",join(",",sort keys %{$info_circs{$id}{reads}}),"\n";
}
close FOUT;

open FOUT,">$outdir/junction_read_details.xls";
print FOUT "read_id\tstrand\tquery1_start\tquery1_end\thit1_start\thit1_end\tquery2_start\tquery2_end\thit2_start\thit2_end\tcircrna_id\tread_type\tNumber of non-encoded nts\tSequence of non-encoded nts\tPercentage of A\tPercentage of T\tPercentage of G\tPercentage of C\n";
foreach my $read_id(sort keys %info_reads){
	foreach my $circ_id(sort keys %{$info_reads{$read_id}{circrna}}){
		if(exists $info_insert{$read_id}{$circ_id}{length}){
			print FOUT $info_reads{$read_id}{circrna}{$circ_id},"\t",$circ_id,"\t",$info_reads{$read_id}{type},"\t",$info_insert{$read_id}{$circ_id}{length},"\t",$info_insert{$read_id}{$circ_id}{sequence},"\t",$info_insert{$read_id}{$circ_id}{percent_A},"\t",$info_insert{$read_id}{$circ_id}{percent_T},"\t",$info_insert{$read_id}{$circ_id}{percent_G},"\t",$info_insert{$read_id}{$circ_id}{percent_C},"\n";
		}else{
			print FOUT $info_reads{$read_id}{circrna}{$circ_id},"\t",$circ_id,"\t",$info_reads{$read_id}{type},"\t\t\t\t\t\t\t\t\t\t\t\t\t\n";
		}
		
	}
}
close FOUT;
