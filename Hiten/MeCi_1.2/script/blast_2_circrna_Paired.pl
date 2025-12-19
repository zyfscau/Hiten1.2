#!/usr/bin/env perl
use strict;
use warnings;

#Author: guantao.zheng
#Create: 20210210

my $indir = $ARGV[0];
my $outdir= $ARGV[1];
my $prefix= $ARGV[2];
my $overlap= $ARGV[3];
my $gap= $ARGV[4];
if(scalar(@ARGV) < 0){
	die "usage: $0 indir outdir !";
}

my $max_len = 30000;

my %convert;
my %info_reads;
my %info_circs;
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
			my $start_a = $t2_s;
			my $start_b = $t2_s + $diff;
			my $end_a   = $t1_e;
			my $end_b   = $t1_e;
			my $length = $end_a - $start_a + 1;
			my $strand = "-";
			
			#for(my $i = $t2_s_a; $i <= $t2_s_b; $i++){
			#	my $x = $i;
			#	my $y = $i + $length - 1;
			#	unless(exists $convert{"$prefix:$x|$y"}){
			#		$convert{"$prefix:$x|$y"} = "$prefix:$t2_s_a|$t2_s|$diff";
			#	}
			#}

			my $id = "$prefix|$end_a|$start_a|$diff|-";
			$info_circs{$id}{start}  = $end_a;
			$info_circs{$id}{end}    = $start_a;
			$info_circs{$id}{strand} = $strand;
			$info_circs{$id}{reads}{$read}++;
			$info_reads{$read}{circrna}{$id} = $line;
		}else{
			next if ($t2_e - $t1_s) > $max_len;
			next if $t2_s <= $t1_e;
			my $start_a = $t2_e;
			my $start_b = $t1_s + $diff;
			my $end_a   = $t1_s;
			my $end_b   = $t2_e;
			my $length = $end_a - $start_a + 1;
			my $strand = "+";

			my $id = "$prefix|$end_a|$start_a|$diff|+";
			$info_circs{$id}{start}  = $end_a;
			$info_circs{$id}{end}    = $start_a;
			$info_circs{$id}{strand} = $strand;
			$info_circs{$id}{reads}{$read}++;
			$info_reads{$read}{circrna}{$id} = $line;
		}
	}
	close FIN;
}

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

open FOUT,">$outdir/reads_details.xls";
print FOUT "read_id\tstrand\tquery1_start\tquery1_end\thit1_start\thit1_end\tquery2_start\tquery2_end\thit2_start\thit2_end\tcircrna_id\tread_type\n";
foreach my $read(sort keys %info_reads){
	foreach my $id(sort keys %{$info_reads{$read}{circrna}}){
		print FOUT $info_reads{$read}{circrna}{$id},"\t",$id,"\t",$info_reads{$read}{type},"\n";
	}
}
close FOUT;
