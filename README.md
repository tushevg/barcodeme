# barcodeme
a script to generate DNA barcodes, given kmer length and Hamming distance

## prerequisite
none

## help
```
barcodeme
version 0.1, Apr 2021

usage:
        $ perl barcodeme.pl [options] > barcodes_list.txt

description:
        generates a list of DNA barcodes, given kmer length and
        Hamming distance, outputs to stdout

options:
-k, --kmer
        size of DNA kmer, restricted [4,10], default 6
-d, --dist
        size of Hamming distance between accepted kmers, restricted [0, kmer], defualt 3
-l, --list
        size of final kmer list, default -1 (all acceptable)
-o, --output
        file name to output a distance matrix of accepted kmers, default skip
-h, --help
        print help message
```

## example
```
perl barcodeme.pl -k 6 -d 3 -l 96 -o barcodes_dist.txt > barcodes_list.txt
```

## references
John H. [Conway](https://ieeexplore.ieee.org/document/1057187) and N. J. A. Sloane<br>
Lexicographic Codes: Error-Correcting Codes from Grame Theory<br>
IEEE Transactions on information theory, Volume: 32, Issue: 3, May 1986<br>

Dan [Ashlock](https://ieeexplore.ieee.org/document/1004430), Ling Guo and Fang Qiu<br>
Greedy Closure Evolutionary Algorithms<br>
IEEE Proceedings of the 2002 Congress on Evolutionary Computation. CEC'02 (Cat. No.02TH8600), August 2002<br>
