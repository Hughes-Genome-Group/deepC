#!usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

# DESCRIPTION
# Script to prune and extract the features of a HiCPro matrix output given a genomic region of interest
# input: 
# 	bin: .bed file with the genomic coords of each bin (id)
# 	matrix: .matrix file (raw or normalized) as output from HiCPro
#	chr:	chromosome of interest (e.g. "chr11")
#	start:	start of genomic region of interst (if not specificied will return the whole chromosome
#	end: end of genomic region of interest
#
# output: pruned .matrix and .bed file (outmatrix and outbed arguments) file; containing only the ids and interaction between the desired coords

# define variables for options
my $chr = "";
my $region_start = 0;
my $region_end = 0;

&GetOptions ( 
	"help|?" => \my $help,
	"bin=s" => \my $bin_coords,
	"matrix=s" => \my $matrix_in,            
	"chr=s" => \$chr,
	"start=i" => \$region_start,
	"end=i" => \$region_end,
	"outbed=s" => \my $out_bed,
	"outmatrix=s" => \my $out_matrix
);

# HELP MESSAGE ---------------------------------------------------------
if($help){
	print "\nScript to extract only interactions within a desired genomic range from HiCPro Matrix output and accorind bin coordinate bed file!\n\n";
	print "Usage: perl extract_hicpro_matrix.pl --bin <bin_coords.bed> --matrix <hicpro.matrix> --outbed <output.coords.bed> --outmatrix <output.pruned.matrix> --chr chrX (--start start_coord_bed --end end_coord_bed)\n\n";
	print "Note: Use 0-based bed coordinates; Will extract the whole chrom if no end coord is supplied!\n\n";
	exit;
}

# declaresome variables
my (@col, %hash);


# 1) Read in bed file, filter for regions selected and store in hash

# check if to retrieve only the chromosome or a efined regions
my $only_chr = 0; 
if($region_end == 0){
	print "No end coordinate selected! Will retrieve the entire chromosome $chr\n";
	$only_chr = 1;
}
print "No start coordinate supplied. WIll retrieve the entire chromosome up to $region_end\n" if ($region_start == 0 && $region_end != 0);


# Start reading in the file
open(OUTBED, ">$out_bed") or die "Can't open output bed file for pruned bin coords.bed $out_bed $1\n";
open(BED, $bin_coords) or die "Can't open Input bin coordinates file $bin_coords $!\n";

	while(<BED>){

		chomp;
		
		@col = split(/\t+/, $_);

		# filter for genomic regions & store in hash
		next if $col[0] ne $chr;

		if($only_chr == 1){	#store ids from entire chrom

				$hash{$col[3]} = 1;  #store id

				print OUTBED $_."\n"; #print column in pruned bed file

				next;

		}
	
		#very simple interset: every bin that starts before or ends after the defined region wil be discarded
		if($col[1] >= $region_start && $col[2] <= $region_end){
		
			$hash{$col[3]} = 1;  #store bin id in hash
		
			print OUTBED $_."\n"; #print column in pruned bed file

			next;

		}

	}

close(BED);
close(OUTBED);


# 2) Read in Matrix file and only print out those entries corresponding to intersactions between stored bin ids
open(OUTMAT, ">$out_matrix") or die "Can't open output matric file for pruned bin matrix $out_matrix $1\n";
open(MAT, $matrix_in) or die "Can't open Input .matrix file $matrix_in $!\n";


	while(<MAT>){

		chomp;

		@col = split(/\t+/, $_);

		#print if start bin and target bin ids both exist in the hash
		if( exists $hash{$col[0]} && exists $hash{$col[1]}){

			print OUTMAT $_."\n";

		}
	}

close(MAT);
close(OUTMAT);

