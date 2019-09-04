use 5.14.0;
use Bio::DB::Fasta;
use Bio::SearchIO;
use Bio::SeqIO;
use File::Copy;
use File::Temp qw(tempfile tempdir);
use Carp;
use Fcntl qw(:flock);

=head1 Name

reciblastn.pl

=head2 BASH script to run this script in parallel instances.

 para_recibl () {
 nj=4;
 nn=$(( nj - 1 ));
 for n in $(seq 0 $nn); do
   echo perl reciblastn.pl $n $nj 900
 done
 }
 rm reciblast.out
 para_recibl | parallel --jobs $nj --dry-run

The 900 above is simply to limit the number of sequences done during testing.
It becomes $testCnt via $ARGV[2] in this script.

=cut


my $tempdir = qw(/tmp);
my $template="jitreciXXXXX";
my $blastbindir = qq(/usr/local/bin);
my $fna1 = qq(sco_prot.fna); # ref. query1
my $fna2 = qq(sven_prot.fna); # sub. query2
my $db = qq(sven_prot); # subject blastdb
my $refdb = qq(sco_prot); # reference blastdb
my $outfile = qq(reciblast.out);
my $biodb = Bio::DB::Fasta->new($fna2);


open(my $ofh, ">>", $outfile);
select($ofh);

my $paranum = $ARGV[0];
my $parajobs = $ARGV[1];
my $testCnt = 0;
if($ARGV[2]) { $testCnt = $ARGV[2]; }


my $seqio = Bio::SeqIO->new(-file => $fna1, -format => 'fasta');
my $seqCnt = 0;
while(my $seqobj = $seqio->next_seq()) {
  $seqCnt += 1;
  if(($seqCnt % $parajobs) == $paranum) {
  my ($forward, $reverse) = reciblastn(query => $seqobj, db => $db,
      refdb => $refdb, biodb => $biodb);
  my $reciprocal = 0;
  if($forward and $reverse) {
    if($forward->{qname} eq $reverse->{hname}) {
      $reciprocal = 1;
    }
    flock($ofh, LOCK_EX);
    print(join("\t",
          $paranum,
          $forward->{qname},
          $forward->{hname},
          $reverse->{hname},
          $forward->{fracid},
          $forward->{signif},
          $reciprocal),
        "\n");
    flock($ofh, LOCK_UN);
  }
  elsif($forward) {
    flock($ofh, LOCK_EX);
    print(join("\t",
          $paranum,
          $forward->{qname},
          $forward->{hname},
          "-",
          $forward->{fracid},
          $forward->{signif},
          $reciprocal),
        "\n");
    flock($ofh, LOCK_UN);
  }
  }
  if($testCnt and $seqCnt >= $testCnt) {
    last;
  }
}

# {{{ sub reciblastn (hash(query, refdb, db, biodb, expect)).
# Returns two hashrefs or one hashref and undef.
# query is a faa file with a single sequence or a protein seqobj.
# db is the subject blast database
# refdb is the reference blast database. i.e. of the organism from which
#       the query comes.
# biodb is the Bio::DB::Fasta object from which the hit id can be retrieved.
# expect is the evalue threshold

sub reciblastn {
  my %args = @_;
  my $query = $args{query};
  my $db = $args{db};
  my $biodb = $args{biodb};
  my $refdb = $args{refdb};
  my $evalue = $args{expect};
  unless ($evalue) { $evalue = 1; }
  my $outfmt = 0;

  my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
  close($fh1);
  if(ref($query)) {
    my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.fna');
    my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    $seqout->write_seq($query);
    close($fh);
    my $xstr = qq($blastbindir/blastn -outfmt $outfmt -query $fn -db $db -evalue $evalue -out $fn1);
    $xstr .= qq( -dust no);
    qx($xstr);
    unlink($fn);
  }
  elsif(-e $query and -r $query) {
    my $xstr = qq($blastbindir/blastn -outfmt $outfmt -query $query -db $db -evalue $evalue -out $fn1);
    $xstr .= qq( -dust no);
    qx($xstr);
  }
# topHSPtopHit
  my %forward = topHSPtopHit($fn1, "blast");
  unlink($fn1);

  my %reverse;
  if(%forward) {
  my $fhname = $forward{hname};
  # carp("Blast.pm: $seqCnt $fhname");
  my $revquery = $biodb->get_Seq_by_id($fhname);
    my($fh2, $fn2)=tempfile($template, DIR => $tempdir, SUFFIX => '.fna');
    my $seqout2 = Bio::SeqIO->new(-fh => $fh2, -format => 'fasta');
    $seqout2->write_seq($revquery);
    close($fh2);
  my($fh3, $fn3)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
    my $qxstr = qq($blastbindir/blastn -outfmt $outfmt -query $fn2 -db $refdb -evalue $evalue -out $fn3);
    $qxstr .= qq( -dust no);
    qx($qxstr);
  %reverse = topHSPtopHit($fn3, "blast");
  unlink($fn2);
  unlink($fn3);
  }
  else { return(); }

  if(%forward and %reverse) {
  return(\%forward, \%reverse);
  }
  elsif(%forward) {
    return(\%forward, undef);
  }
  else {
    return();
  }
}
# }}}


# {{{ topHSPtopHit (blastOutputFileName) returns(list of hashes(qname, hname, qlen, hlen, signif, bit hdesc, qcover,
# hcover, hstrand) );
# Gives the top HSP of only the top hit in each blast result in a file.
sub topHSPtopHit {
  my $filename=shift(@_);
  my $format = 'blast';
  my $temp = shift(@_);
  if($temp) { $format = $temp; }
#print(STDERR "in topHit $filename\n");
  my $searchio = Bio::SearchIO->new(-format => $format,
      -file => $filename
      );
  my %rethash;
  my $result = $searchio->next_result();
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qlen=$result->query_length();
  my $hit = $result->next_hit();
  if ($hit) {
    my $num_hsps = $hit->num_hsps();
    my $hsp = $hit->next_hsp();
    if ($hsp) {
      my $hname=$hit->name();
      my $hlen=$hit->length();
      my $frac_id = sprintf("%.3f", $hsp->frac_identical());
      my $hdesc=$hit->description();
      my $signif=$hsp->significance();
      my $laq=$hsp->length('query');
      my $lah=$hsp->length('hit');
      my $qcov = sprintf("%.3f", $laq/$qlen);
      my $hcov = sprintf("%.3f", $lah/$hlen);
      my $qstart = $hsp->start('query');
      my $qgaps = $hsp->gaps("query");
      my $hgaps = $hsp->gaps("hit");
      my $qend = $hsp->end('query');
      my $hstart = $hsp->start('hit');
      my $hend = $hsp->end('hit');
      my $hframe = $hsp->frame('hit');
      my $bitScore = $hsp->bits();
      my $strand = $hsp->strand('hit');
      %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
          signif => $signif, bit => $bitScore, hdesc => $hdesc,
          hstrand => $strand, qstart => $qstart, hframe => $hframe,
          qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
          fracid => $frac_id, qcover => $qcov, qcov => $qcov, hcov => $hcov,
          hcover => $hcov, numhsps => $num_hsps, qgaps => $qgaps, hgaps => $hgaps);
    }
  }
  return(%rethash);
}
# }}}

