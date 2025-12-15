#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Config::IniFiles;
use Cwd qw/abs_path/;
use FindBin qw($Bin $Script $RealBin);
use File::Basename qw(basename dirname);
use Parallel::ForkManager;
use POSIX;

my $script = "$Bin/script";
my $bin    = "$Bin/bin";

my ($fq_1, $fq_2, $genome, $chunksize, $outdir, $help);
my ($max_process, $blast_cpu, $blast_evalue);
my $Paired;
my $output_prefix;
my $overlap;
my $gap;

GetOptions(
	"in1:s"           => \$fq_1,
	"in2:s"           => \$fq_2,
	"genome:s"        => \$genome,
	"chunksize:i"     => \$chunksize,

	"outdir:s"        => \$outdir,

	"max_process:i"   => \$max_process,
	"blast_cpu:i"     => \$blast_cpu,
	"blast_evalue:s"  => \$blast_evalue,
	"output_prefix:s" => \$output_prefix,
	"overlap:i"       => \$overlap,
	"gap:i"           => \$gap,	
	
	"help!"           => \$help,
	"Paired"          => \$Paired,	
);

&parsing_parameters();

my @fa_files;
my %IDs;
my $out_index = "$outdir/01.index";
my $out_reads = "$outdir/02.reads";
my $out_split = "$outdir/03.split";
my $out_blast = "$outdir/04.blast";
my $out_parse = "$outdir/05.parse";


system("mkdir -p $outdir $out_index $out_reads $out_split $out_blast $out_parse");

&index_genome();
&parse_reads();
&run_split();
&run_blast_and_parse();

print "All jobs done!\n";

###########################################################################

sub parsing_parameters{
	if($help){
		&usage;
	}
	if(!defined($fq_1)){
		print STDERR "please specifiy in1!";
		&usage;
	}elsif(! -e $fq_1){
		print STDERR "in1 file[$fq_1] does not exist!";
		&usage;
	}
	if(defined($fq_2)){
		if(! -e $fq_2){
			print STDERR "in2 file[$fq_2] does not exist!";
			&usage;
		}
	}
	
	unless(defined($genome)){
		print STDERR "please specifiy genome file!";
		&usage;
	}
	
	$chunksize      ||= 1000000;
	$blast_cpu      ||= 10;
	$blast_evalue   ||= 2;
	$max_process    ||= 4;
	$overlap        ||= 3;
	$gap            ||= 3;
	$outdir         ||= "out";
	$output_prefix  ||= "PL";
	
	if($chunksize <= 0){
		die "Error: chunksize must larger than 0 !!\n";
	}
	
	system("mkdir -p $outdir");
	$outdir       = abs_path($outdir);
	$genome       = abs_path($genome);
	$fq_1         = abs_path($fq_1);
	if($fq_2){
		$fq_2     = abs_path($fq_2);
	}
}
sub parse_reads{
	print STDERR "run parse reads...\n";
	if(defined($fq_2)){
		# read 2 : sense 
		system("$bin/flash $fq_2 $fq_1 --threads 10 --output-prefix merged --max-overlap 150 --output-directory $out_reads &> $out_reads/merged.log");
		#system("ln -sf $out_reads/merged.extendedFrags.fastq $out_reads/merged.fq");
		# system("$bin/fastqToFa $out_reads/merged.fq $out_reads/reads.fa");
		system("ln -sf $out_reads/merged.extendedFrags.fastq $out_reads/part_merged.fq");
		system("$script/modify_fastq.pl -i $out_reads/merged.notCombined_1.fastq -o $out_reads/part_2.fq --read 2");
		system("$script/modify_fastq.pl -i $out_reads/merged.notCombined_2.fastq -o $out_reads/part_1.fq --reverse --read 1");
		system("cat $out_reads/part_merged.fq $out_reads/part_1.fq $out_reads/part_2.fq > $out_reads/reads.fq");
		system("$bin/fastqToFa $out_reads/reads.fq $out_reads/reads.fa");
	}else{
		system("$bin/fastqToFa $fq_1 $out_reads/reads.fa");
	}
}

sub index_genome{
	print STDERR "run index genome...\n";
	system("ln -sf $genome $out_index/genome.fa");
	system("$bin/formatdb -i $out_index/genome.fa -n $out_index/index_blast -t index_blast -l $out_index/index_blast.log -p F");
}

sub run_split{
	print STDERR "run fasta split...\n";
	my $seq_file = "$out_reads/reads.fa";
	my $seq_num  = fastalength($seq_file);
	my $file_num = POSIX::ceil($seq_num/$chunksize);
	for(my $i = 1; $i <= $file_num; $i++){
		my $x = convert_id($i);
		$IDs{convert_id($i)}++;
		push @fa_files, "$out_split/reads.fa_$x.split";
	}
	
	my $fa_file = "$out_split/reads.fa";
	
	system("ln -sf $seq_file $fa_file");
	system("cd $out_split/; $script/splitfasta.pl $fa_file $chunksize; cd -" );
}

sub fastalength {
	my $fasta = shift;
	my $num   = `grep '^>' -c $fasta`;
	chomp($num);

	return $num;
}

sub run_blast_and_parse {
	print STDERR "run run_blast_and_parse...\n";
	my $pm = new Parallel::ForkManager($max_process);
	foreach my $fa_file (@fa_files) {
		my $pid = $pm->start and next;
		my $fa_name    = basename($fa_file);
		
		system("$bin/fastalength $fa_file \| awk -F \' \' \'{print \$2\"\\t\"\$1}\' \> $out_split/$fa_name.len");
		system("$bin/blastall -p blastn -d $out_index/index_blast -i $out_split/$fa_name -a $blast_cpu -e $blast_evalue -m 8 -o $out_blast/$fa_name.txt -F F \> $out_blast/$fa_name.log");
		system("$script/parsing_blast.pl $out_split/$fa_name.len $out_blast/$fa_name.txt $out_parse/$fa_name.txt");
		
		$pm->finish;
	}
	$pm->wait_all_children;
	if(defined($Paired)){
		system("$script/blast_2_circrna_Paired.pl $out_parse $outdir $output_prefix $overlap $gap");
		print STDERR "Paired end\n";
	}else{
		system("$script/blast_2_circrna_Single.pl $out_parse $outdir $output_prefix $overlap $gap $out_reads/reads.fa");
		print STDERR "Single end\n";
	}
}

sub convert_id{
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

sub usage{
	my $program = basename($0);
    die "
Version: v1.2
Create : 20220301
Update : 20240303
Aurthor: guantao.zheng\@origin-gene.com
Usage  : $program <arguments>
Arguments:
  Input Options:
    --in1           <string>    read1 input file name in fastq format, [required]
    --in2           <string>    read2 input file name in fastq format, [required for PE input]
    --genome        <string>    genome file in fasta format, [ required ]
    --chunksize     <number>    chunk size for sequence split, default is 1000000 [ optional ]

  Advanced Options:
    --blast_cpu     <number>    blast cpu number. default is 10  [ optional ]
    --blast_evalue  <string>    blast evalue threshold. default is 2  [ optional ]
    --max_process   <number>    max process number limitation for sequence alignment, default is 4 [ optional ]
    --Paired        <null>      Paired-end mode. the input file is paired-end [ optional ]
    --output_prefix <string>    Prefix of output circRNA ID.  Default is 'PL' [ optional ]
  General Options:
    --outdir        <string>    output dir ,default is 'out'
    --help          <null>      display help information [ optional ]
\n";
}
