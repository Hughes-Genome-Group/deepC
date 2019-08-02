# perl script for matching a hic interaction matrix to a query table
# query and hic coord left most must be the smallest

my $hic_in = $ARGV[0];
my $query_in = $ARGV[1];
my $out = $ARGV[2];
my $binsize = $ARGV[3];

my (%h, $l, $r);
my $zero = 0;

# 1) Read HiC data into greedy perl hash
open(HI, $hic_in) or die "Can't open HiC input matrix $hic_in $!\n";

  while(<HI>){

    chomp;
    my ($x, $y, $v) = split(/\t+/, $_);

    # get leftmost and rightmost coordinate
    if($x >= $y){
      $l = $y;
      $r = $x;
    }else{
      $l = $x;
      $r = $y;
    }
    # adjust to leftmost bin notation ### IS THAT CORRECT?!
    $l -= $binsize/2;
    $r -= $binsize/2;

    my $tag = "$l:$r";

    if(! exists $h{$tag}){
      $h{$tag} = $v;
    }else{
      print("Duplicate Value! $tag ...\n");
    }

  }

close(HI);



# 2) Run through query table --> grep queired interactions and print
open(IN, $query_in) or die "Can't open query table input $query_in $!\n";
open(OUT, ">$out") or die "Can't open output file $out $!\n";

  while(<IN>){

    chomp;
    my ($chr, $start, $end, @queries) = split(/\t+/, $_);

    # foreach query tag --> query the interaction value from the hash
    # set 0 if no interaction found
    my @values;

    foreach my $q (@queries){
      if(exists $h{"$q"}){
        push(@values, $h{"$q"});
      }else{
        push(@values, $zero);
      }
    }

    # print line
    print OUT "$chr\t$start\t$end\t".join("\t", @values)."\n";

  }

close(IN);
close(OUT);

print "Finished ...\n";
