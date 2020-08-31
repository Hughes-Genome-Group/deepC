# perl script for creating the query table for ma zigzag vertical pole hic encoding

my $binned_in = $ARGV[0];
my $binsize = $ARGV[1];

my $outfile = $ARGV[2];

my $halfbin = $binsize/2;

# Read binned genome and make zigzag query pole over centerh
open(IN, $binned_in) or die "Can't open binned genome file $binned_in $!\n";
open(OUT, ">$outfile") or die "Can't open output file $outfile $!\n";
	while(<IN>){

	chomp;
	my ($chr, $start, $end) = split(/\t+/, $_);

    	# get center position
    	my $center = ($end - $start)/2 + $start;
	my $center_left_adjust = $center - $halfbin;
	my $end_left_adjust = $end - $bin;

	# lay out queries (x)
	my @x;
	my @y;

	for(my $i=$start; $i <= $center_left_adjust; $i = $i + $binsize){
		# repeat every element twice
		push(@x, $i);
		push(@x, $i);
	}

	#print join("\t", @x)."\n";
        for(my $i=$center_left_adjust; $i < $end_left_adjust; $i = $i + $binsize){
                # repeat every element twice
                push(@y, $i);
                push(@y, $i);
        }

	# trim last positions
	pop(@x);
	pop(@y);

	# revert y
	@y = reverse @y;

	# get x length
	my $x_length = @x;

	# combine to queries
	my @q;
	for(my $j=0; $j < $x_length; $j++){

		my $tmp = $x[$j].":".$y[$j];
		push(@q, $tmp);

	}

	@q = reverse @q;

	# print
	print OUT "$chr\t$start\t$end\t".join("\t", @q)."\n";

	}

close(IN);
close(OUT);
