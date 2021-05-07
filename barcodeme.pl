#!/usr/bin/perl

# barcodeme
# a script to generate DNA barcodes
# given kmer length and Hamming distance
#
# author: Georgi Tushev
# Scientific Computing Facility
# Max Planck Institute for Brain Research
# bugs: sciclist@brain.mpg.de
# v0.1, Apr 2021

use warnings;
use strict;
use Getopt::Long;

sub usage($);
sub parseArgs();
sub isRepeat($);
sub isPalindrome($);
sub isGCbiased($$$);
sub editDistance($$);
sub indexIterator($$);
sub indexToKmer($$);
sub generatePool($$$);
sub findSet($$$);
sub printSet($$);
sub printDist($$);

main:{
    my ($kmer_size,
        $dist_value,
        $list_size,
        $file_table) = parseArgs();

    my @alphabet = qw/A C G T/;
    my @kmer_pool = ();
    my @kmer_set = ();
    
    generatePool(\@kmer_pool, $kmer_size, \@alphabet);
    print STDERR "pool size: ", scalar(@kmer_pool), "\n";

    # magic seed word is the zero word
    my $magic_seed = ($alphabet[0]) x $kmer_size;
    push(@kmer_set, $magic_seed);

    findSet(\@kmer_set, \@kmer_pool, $dist_value);
    shift @kmer_set; # remove magic seed
    print STDERR "set size: ", scalar(@kmer_set), "\n";

    printSet(\@kmer_set, $list_size);
    printDist(\@kmer_set, $file_table) if (defined($file_table));
}



## printDist
## prints distance matrix to file
sub printDist($$)
{
    my $kmer_list = $_[0];
    my $file_table = $_[1];

    open(my $fh, ">", $file_table) or die $!;

    foreach my $kmer_qry (@{$kmer_list}) {
        my @dist_list = ();
        foreach my $kmer_ref (@{$kmer_list}) {
            my $dist = editDistance($kmer_qry, $kmer_ref);
            push(@dist_list, $dist);
        }
        print $fh $kmer_qry, "\t", join(",", @dist_list), "\n";
    }

    close($fh);
}



## printSet
## prints final set to STDOUT
sub printSet($$)
{
    my $kmer_set = $_[0];
    my $list_size = $_[1];

    if ((0 < $list_size) && ($list_size <= scalar(@{$kmer_set}))) {
        @{$kmer_set} = @{$kmer_set}[0 .. ($list_size - 1)];
    }

    print join("\n", @{$kmer_set}), "\n";
}



## findSet
## Conway’s Lexicode Algorithm
sub findSet($$$)
{
    my $kmer_set = $_[0];
    my $kmer_pool = $_[1];
    my $dist_value = $_[2];
    foreach my $kmer_qry (@{$kmer_pool}) {
        my $dist_max = length($kmer_qry);
        foreach my $kmer_ref (@{$kmer_set}) {
            my $dist = editDistance($kmer_qry, $kmer_ref);
            $dist_max = $dist if ($dist_max >= $dist);
        }
        push(@{$kmer_set}, $kmer_qry) if ($dist_value <= $dist_max);
    }
}



## generatePool
## explore all permutation of given alphabet
## filter repeats is crucial for the magic seed ;)
sub generatePool($$$)
{
    my $kmer_pool = $_[0];
    my $kmer_size = $_[1];
    my $alphabet = $_[2];
    my $iterator = indexIterator($kmer_size, scalar(@{$alphabet}));
    while (my $index = $iterator->()) {
        my $kmer = indexToKmer($index, $alphabet);
        next if (isRepeat($kmer));
        next if (isPalindrome($kmer));
        next if (isGCbiased($kmer, 0.3, 0.7));
        push(@{$kmer_pool}, $kmer);
    }
}



## indexToKmer
## convert base 4 index to DNA
sub indexToKmer($$)
{
    my $index = $_[0];
    my $alphabet = $_[1];
    return join("", map { $alphabet->[$_] } @{$index});
}



## indexIterator
## all permutation for given
## word size and alphabet size
sub indexIterator($$)
{
    my $k = $_[0];
    my $alphabet_size = $_[1];
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
                last if ($index[$i] < $alphabet_size);
                $index[$i] = 0;
                $step-- if ($step == $i);
            }
        }
        return undef if ($step < 0);
        return \@index;
    };
}



## editDistance
## Hamming distance (number of mismatches)
sub editDistance($$)
{
    my $dist = $_[0] ^ $_[1];
    return ($dist =~ tr/\0//c);
}



## isRepeat
## use RegEx to filter various seq. repeats 
sub isRepeat($)
{
    my $kmer = $_[0];
    my $is_repeat = 0;

    # skip mono-repeats (> 2)
    $is_repeat = 1 if ($kmer =~ m/([ACGT])\1{2,}/); 

    # skip 2-repeats (e.g. ACACACAC) 
    $is_repeat = 1 if ($kmer =~ m/([ACGT][ACGT])\1{2,}/);

    # skip 3-repeats (e.g. ACGACGACG)
    $is_repeat = 1 if ($kmer =~ m/([ACGT][ACGT][ACGT])\1{1,}/);

    return $is_repeat;
}



## isPalindrome
## forward equals backwards
sub isPalindrome($)
{
    return ($_[0] eq reverse($_[0]));
}



## isGCbiased
## accepted range [0.3, 0.7]
## can be added as user argument
sub isGCbiased($$$)
{
    my $kmer = $_[0];
    my $gc_min = $_[1];
    my $gc_max = $_[2];
    my $gc = ($kmer =~ tr/[GC]//c) / length($kmer);
    return (($gc < $gc_min) || ($gc_max < $gc));
}



## parseArgs
## parses @ARGV for options
sub parseArgs()
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
