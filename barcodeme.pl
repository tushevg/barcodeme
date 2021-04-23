#!/usr/bin/perl

# barcodeme
# a script to generate DNA barcodes
# given kmer length and Hamming distance
#
# author: Georgi Tushev
# bugs: sciclist@brain.mpg.de
# v0.1, Apr 2021

use warnings;
use strict;
use Getopt::Long;

sub usage($);
sub parseargs();
sub kmersampler($);
sub samplepool($$);
sub hammingdist($$);
sub findset($$$);
sub printset($$$);
sub printdist($$);

main:{
    my ($kmer_size,
        $dist_value,
        $list_size,
        $file_table) = parseargs();

    my %kmer_pool = ();
    my %kmer_set = ();
    my @kmer_list = ();

    samplepool(\%kmer_pool, $kmer_size);
    findset(\%kmer_set, \%kmer_pool, $dist_value);
    printset(\@kmer_list, \%kmer_set, $list_size);
    printdist(\@kmer_list, $file_table) if (defined($file_table));
}



## printdist
## prints distance matrix to file
sub printdist($$)
{
    my $kmer_list = $_[0];
    my $file_table = $_[1];

    open(my $fh, ">", $file_table) or die $!;

    foreach my $kmer_qry (@{$kmer_list}) {
        my @dist_list = ();
        foreach my $kmer_ref (@{$kmer_list}) {
            my $dist = hammingdist($kmer_qry, $kmer_ref);
            push(@dist_list, $dist);
        }
        print $fh $kmer_qry, "\t", join(",", @dist_list), "\n";
    }

    close($fh);
}



## printset
## prints final set to STDOUT
sub printset($$$)
{
    my $kmer_list = $_[0];
    my $kmer_set = $_[1];
    my $list_size = $_[2];

    foreach my $kmer (sort keys %{$kmer_set}) {
        print $kmer, "\n";
        push(@{$kmer_list}, $kmer);
        $list_size--;
        last if ($list_size == 0);
    }
}



## findset
## randomly seeds a hash set
## extends set with given distance
sub findset($$$)
{
    my $kmer_set = $_[0];
    my $kmer_pool = $_[1];
    my $dist_value = $_[2];

    # key step is choosing seed kmer
    # greedy algorithm for random seed in pool
    # everytime generates a different set
    print STDERR "find kmer set ... ";

    foreach my $kmer_qry (keys %{$kmer_pool}) {
        my $dist_min = length($kmer_qry);
        foreach my $kmer_ref (keys %{$kmer_set}) {
            my $dist = hammingdist($kmer_qry, $kmer_ref);
            $dist_min = $dist if($dist < $dist_min);
        }

        $kmer_set->{$kmer_qry}++ if($dist_min >= $dist_value);
    }

    print STDERR "done. size: ", scalar(keys %{$kmer_set}), "\n";
}



## hammingdist
## count mismatches in strings with equal length
sub hammingdist($$)
{
    my $word_a = $_[0];
    my $word_b = $_[1];

    return -1 if(length($word_a) != length($word_b));

    my $dist = $word_a ^ $word_b;
    return $dist =~ tr/\0//c;
}



## samplepool
## iterate all kmer permutations
sub samplepool($$)
{
    my $kmer_pool = $_[0];
    my $kmer_size = $_[1];
    my $iterator = kmersampler($kmer_size);
    
    print STDERR "sample kmer pool ... ";

    while (my $kmer = $iterator->()) {
    
        # skip kmers with repeated bases (> 2)
        next if($kmer =~ m/([ACGT])\1{2,}/); 
        
        # skip palindroms
        next if ($kmer eq reverse($kmer));
        
        # add more filters (e.g. GC-content)
        $kmer_pool->{$kmer}++;
    }
    
    print STDERR "done. size: ", scalar(keys %{$kmer_pool}), "\n";
}



## kmersampler
## iterator over kmer permutations
sub kmersampler($)
{
    my $k = $_[0];
    my @alphabet = ("A", "C", "T", "G");  
    my $alphabet_size = scalar(@alphabet);
    my @index;
    my $step = ($k - 1);

    return sub {
        if (scalar(@index) == 0) {
            @index = (0) x $k;
        }
        else {
            my $i = $k;
            while (--$i >= 0) {
                $index[$i]++;
                last if($index[$i] < $alphabet_size);
                $index[$i] = 0;
                $step-- if ($step == $i);
            }
        }
        return undef if($step < 0);
        return join("", map { $alphabet[$_] } @index);
    };
}



## parseargs
## parses @ARGV for options
sub parseargs()
{
    # defaults
    my $kmer_size = 6;
    my $dist_value = 3;
    my $list_size = -1;
    my $file_table;
    my $help;

    Getopt::Long::GetOptions(
        "kmer|k=i" => \$kmer_size,
        "dist|d=i" => \$dist_value,
        "list|l=i" => \$list_size,
        "output|o=s" => \$file_table,
        "help|h" => \$help
    ) or usage("error, unknown command line option");

    if ($help) {
        usage("barcodeme\nversion 0.1, Apr 2021");
    }
    
    if (($dist_value < 0) || ($kmer_size < $dist_value)) {
        usage("error, distance should be non-negative number smaller than kmer size");
    }

    if (($kmer_size < 4) || (10 < $kmer_size)) {
        usage("error, kmer size is restricted [4, 10], devault value 6");
    }
    
    return ($kmer_size, $dist_value, $list_size, $file_table);
}



## usage
## return help message
sub usage($)
{
    my $message = $_[0];
    if (defined($message) && length($message)) {
        $message .= "\n" unless($message =~ /\n$/);
    }
    
    my $command = $0;
    $command =~ s#^.*/##;
    
    print STDERR (
    $message,"\n",
    "usage:\n" .
    "\t\$ perl $command [options] > barcodes_list.txt\n\n" .
    "description:\n" .
    "\tgenerates a list of DNA barcodes, given kmer length and\n" .
    "\tHamming distance, outputs to stdout\n\n" .
    "options:\n" .
    "-k, --kmer\n" .
    "\tsize of DNA kmer, restricted [4,10], default 6\n" .
    "-d, --dist\n" .
    "\tsize of Hamming distance between accepted kmers, restricted [0, kmer], defualt 3\n" .
    "-l, --list\n" .
    "\tsize of final kmer list, default -1 (all acceptable)\n" .
    "-o, --output\n" .
    "\tfile name to output a distance matrix of accepted kmers, default skip\n" .
    "-h, --help\n" .
    "\tprint help message\n"
    );
    
    die("\n");
}
