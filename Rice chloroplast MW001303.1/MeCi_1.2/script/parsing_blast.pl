#!/usr/bin/env perl
use strict;
use warnings;

#Author: guantao.zheng
#Create: 20210210
#
my $file_len   = $ARGV[0];
my $file_blast = $ARGV[1];
my $file_parse = $ARGV[2];

if(scalar @ARGV < 3){
	die "usage: $0 file_len file_blast file_out !";
}

my %length;
open  FIN,"$file_len" or die "can not open $file_len!";
while(my $line = <FIN>){
	chomp($line);
	my @p = split /\t/,$line;
	$length{$p[0]} = $p[1];
}
close FIN;

my %black_list;
my $name;
my %info;
my %ssss;
open FIN,"$file_blast" or die "can not open $file_blast!";
while( my $line = <FIN>){
	my ($query,$chr,$identity,$len,$mis,$gap,$query_s,$query_e,$start,$end,$evalue,$score) = split /\t/,$line;

	if( $length{$query} == $len){
		$black_list{$query}++;
	}
	next if exists $black_list{$query};

	$ssss{$query}++;
	my $num = $ssss{$query};
	my $strand = "";
	if($start < $end){
		$strand = "+";
	}else{
		$strand = "-";
		my $tmp = $end;
		$end = $start;
		$start = $tmp;
	}
	$info{$query}{$num}{strand} = $strand;
	$info{$query}{$num}{q_s}    = $query_s;
	$info{$query}{$num}{q_e}    = $query_e;
	$info{$query}{$num}{t_s}    = $start;
	$info{$query}{$num}{t_e}    = $end;
}
close FIN;

open FOUT,">$file_parse" or die "can not open $file_parse!";
foreach my $query(sort keys %info){
	my @s_p;
	my @s_m;
	foreach my $num( sort { $a<=>$b } keys %{$info{$query}} ){
		if($info{$query}{$num}{strand} eq "+"){
			push @s_p, $num;
		}else{
			push @s_m, $num;
		}
	}
	
	if(scalar(@s_p) >= 2){
		for(my $i=0; $i<(scalar(@s_p)-1); $i++){
			foreach(my $j = $i + 1; $j < scalar(@s_p); $j++){
				if(    $info{$query}{$s_p[$i]}{q_s} == 1 && $info{$query}{$s_p[$j]}{q_e} == $length{$query} && $info{$query}{$s_p[$i]}{t_s} > $info{$query}{$s_p[$j]}{t_s}){
					print FOUT $query,"\t";
					print FOUT "+\t";
					print FOUT $info{$query}{$s_p[$i]}{q_s},"\t";
					print FOUT $info{$query}{$s_p[$i]}{q_e},"\t";
					print FOUT $info{$query}{$s_p[$i]}{t_s},"\t";
					print FOUT $info{$query}{$s_p[$i]}{t_e},"\t";
					print FOUT $info{$query}{$s_p[$j]}{q_s},"\t";
					print FOUT $info{$query}{$s_p[$j]}{q_e},"\t";
					print FOUT $info{$query}{$s_p[$j]}{t_s},"\t";
					print FOUT $info{$query}{$s_p[$j]}{t_e},"\n";
				}elsif($info{$query}{$s_p[$j]}{q_s} == 1 && $info{$query}{$s_p[$i]}{q_e} == $length{$query} && $info{$query}{$s_p[$j]}{t_s} > $info{$query}{$s_p[$i]}{t_s}){
					print FOUT $query,"\t";
					print FOUT "+\t";
					print FOUT $info{$query}{$s_p[$j]}{q_s},"\t";
					print FOUT $info{$query}{$s_p[$j]}{q_e},"\t";
					print FOUT $info{$query}{$s_p[$j]}{t_s},"\t";
					print FOUT $info{$query}{$s_p[$j]}{t_e},"\t";
					print FOUT $info{$query}{$s_p[$i]}{q_s},"\t";
					print FOUT $info{$query}{$s_p[$i]}{q_e},"\t";
					print FOUT $info{$query}{$s_p[$i]}{t_s},"\t";
					print FOUT $info{$query}{$s_p[$i]}{t_e},"\n";
				}
			}
		}
	}
	if(scalar(@s_m) >= 2){
		for(my $i=0; $i<(scalar(@s_m)-1); $i++){
			foreach(my $j = $i + 1; $j < scalar(@s_m); $j++){
				if(    $info{$query}{$s_m[$i]}{q_s} == 1 && $info{$query}{$s_m[$j]}{q_e} == $length{$query} && $info{$query}{$s_m[$i]}{t_s} < $info{$query}{$s_m[$j]}{t_s}){
					print FOUT $query,"\t";
					print FOUT "-\t";
					print FOUT $info{$query}{$s_m[$i]}{q_s},"\t";
					print FOUT $info{$query}{$s_m[$i]}{q_e},"\t";
					print FOUT $info{$query}{$s_m[$i]}{t_s},"\t";
					print FOUT $info{$query}{$s_m[$i]}{t_e},"\t";
					print FOUT $info{$query}{$s_m[$j]}{q_s},"\t";
					print FOUT $info{$query}{$s_m[$j]}{q_e},"\t";
					print FOUT $info{$query}{$s_m[$j]}{t_s},"\t";
					print FOUT $info{$query}{$s_m[$j]}{t_e},"\n";
				}elsif($info{$query}{$s_m[$j]}{q_s} == 1 && $info{$query}{$s_m[$i]}{q_e} == $length{$query} && $info{$query}{$s_m[$j]}{t_s} < $info{$query}{$s_m[$i]}{t_s}){
					print FOUT $query,"\t";
					print FOUT "-\t";
					print FOUT $info{$query}{$s_m[$j]}{q_s},"\t";
					print FOUT $info{$query}{$s_m[$j]}{q_e},"\t";
					print FOUT $info{$query}{$s_m[$j]}{t_s},"\t";
					print FOUT $info{$query}{$s_m[$j]}{t_e},"\t";
					print FOUT $info{$query}{$s_m[$i]}{q_s},"\t";
					print FOUT $info{$query}{$s_m[$i]}{q_e},"\t";
					print FOUT $info{$query}{$s_m[$i]}{t_s},"\t";
					print FOUT $info{$query}{$s_m[$i]}{t_e},"\n";
				}
			}
		}
	}
}
close FOUT;
