##### getUmi.pl 
##### Julie Zhu
##### Oct 26th 2016

#!usr/bin/perl
use strict;

if ($#ARGV <2) {
	print "usage: perl getUmi.pl sequenceDir barcodeLength UMIlength\n";
	print "example: perl getUmi.pl fastq 16 8\n";
	exit(0);
}

my $seqdir = $ARGV[0];
my $barcodeLen = $ARGV[1];
my $UMIlen = $ARGV[2];
my @seqfiles = <$seqdir/*>;
foreach my $file (sort@seqfiles) {
       getUmi($file);
}

sub getUmi 
{
	my $seqfile = shift;
	open(INFILE, "$seqfile") || die("Cannot open file $seqfile  $!\n");
	my @temp = split(/\//,$seqfile);	
	my $outfile = $temp[length(@temp)];
	my $outfile1 = "UMI-" . $outfile;
	my ($line, @line, $umi, $firstbps, $seqlen);
	open(OUTFILE1, ">$outfile1") || die("Cannot open file $outfile1 $!\n");
	my $i =0;
	while (<INFILE>)
	{
		$line = $_;
		chomp $line;
		$i++;
		if ($i % 4 ==1)
		{
		      @line = split(" ", $line);	
		      @temp = split(":", $line[1]);
		      ### extract umi
		      $umi = substr($temp[scalar @temp - 1],$barcodeLen, $UMIlen);
		      print OUTFILE1 "$line[0]\t$umi\n";
		}
	}
	close(INFILE);
	close(OUTFILE1);
}
