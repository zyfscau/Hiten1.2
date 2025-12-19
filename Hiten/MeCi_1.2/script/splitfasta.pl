#!/usr/bin/env perl
use strict;
use warnings;
use Carp;
use File::Basename;
#print basename($0)."\n";
if (@ARGV!=2) {
    die "Usage:splitfasta.pl inputfile fastaseqnumber\n";
}
open(SEQ,"< $ARGV[0]")||die "Can't open $ARGV[0]\n" ;
my $a=0;
my $i=1;
my $file=basename($ARGV[0]);
open(WF,"> $file\_".convert($i).".split")||die"Can't write to $file\_1:$!\n";
while (<SEQ>) {
  chomp;
  next if /^\s*$/;
  if(/^>\s*(\S+)/){
  	my $name=$1;
  	 if($a>=$ARGV[1]){
  	   close WF;
  	   $a=0;
  	   $i++;
  	   open(WF,"> $file\_".convert($i).".split")||die"Can't write to $file\_$i:$!\n";
     }
     $a++;
  	#print "$a\n";
  	print WF ">$name\n";
  }else{
  	print WF "$_\n";
  }
}
close WF;
close SEQ;

sub convert{
	my ($x) = @_;
	if($x<10){
		return "00$x";
	}elsif($x<100){
		return "0$x";
	}elsif($x<1000){
		return "$x";
	}else{
		die "the number of splited files is too much!!! please modity chunk size!!";
	}
}
