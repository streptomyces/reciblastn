# reciblastn.pl


## BASH code to run this reciblastn.pl in parallel instances.

~~~ {.sh}
 para_recibl () {
 nj=4;
 nn=$(( nj - 1 ));
 for n in $(seq 0 $nn); do
   echo perl reciblastn.pl $n $nj 900
 done
 }
 rm reciblast.out
 para_recibl | parallel --jobs $nj --dry-run
~~~

The 900 above is simply to limit the number of sequences done during testing.
It becomes `$testCnt` via `$ARGV[2]` in this script.

## Running

GNU Parallel is needed to run multiple instances of the Perl script. It should
be installed and working.

This is a proof of concept script. There is a lot of hardcoding inside. Sven and
Sco data files are provided for the sake of testing.

## Author

govind.chandra@jic.ac.uk

