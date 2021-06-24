# Make a SNP info file for MPH from a file of genomic feature

# All postions are 1-based.
# start and end are included in invervals.

use strict;
use warnings;

@ARGV == 3 or &print_usage();
my ($snp_pos_file, $gf_file, $output_prefix) = @ARGV;

# step 1
# get SNP positions
open IN,$snp_pos_file or die "Cannot open file $snp_pos_file: $!\n";
my %snp;
my @snp_arr;
while(<IN>)
{
	chomp;
	my @c = split /\s+/;
	push @snp_arr, $c[1];
	push @{$snp{$c[0]}},[@c[2,1]];
}
close IN;

for (keys %snp)
{
	@{$snp{$_}} = sort {$a->[0] <=> $b->[0]} @{$snp{$_}};
}

# step 2
# get genomic feature
open IN,$gf_file or die "Cannot open file $gf_file: $!\n";
my %gf;
while(<IN>)
{
	chomp;
	my @c = split /\s+/;
	push @{$gf{$c[0]}},[@c[1..3]];
}
close IN;

for (keys %gf)
{
	@{$gf{$_}} = sort {$a->[0] <=> $b->[0]} @{$gf{$_}};
}

# step 3
# assign SNPs to genomic features with single pass
my %snp2gf;
for my $chrom (keys %snp)
{
	next unless (defined $gf{$chrom});
	my @snp_chr_loc = @{$snp{$chrom}};
	my @gf_chr_loc = @{$gf{$chrom}};
  
	for (@snp_chr_loc)
	{
		my ($snp_pos,$snp_id) = @$_;
		
		for (my $i=0;$i<@gf_chr_loc;$i++)
		{
			my ($gf_s, $gf_e, $gf_c) = @{$gf_chr_loc[$i]};
			
			if($gf_e < $snp_pos)
			{
				splice @gf_chr_loc,$i,1;
				$i --;
			}
			elsif($gf_s > $snp_pos)
			{
				last;
			}
			else
			{
				$snp2gf{$snp_id} = $gf_c;
				last;
			}
		}
	}
}

# step 4
# generate output file
my @categories = &uniq(values %snp2gf);
my %cat2line;
for (0..$#categories) {
	my @cc = ("") x scalar(@categories);
	$cc[$_] = 1;
	$cat2line{$categories[$_]} = join ",", @cc;
}
my $null_str = "," x $#categories;
open OUT, ">$output_prefix.snp_info.csv";
print OUT join(",", "SNP", "intercept", @categories), "\n";
for (@snp_arr) {
	if(defined $snp2gf{$_}) {
		print OUT join(",", $_, 1, $cat2line{$snp2gf{$_}}), "\n";
	} else {
		print OUT join(",", $_, 1, $null_str), "\n";
	}
}

sub print_usage()
{
	print "Run the program:\n";
	print "  perl enrichment.pl plink-bim-file genomic-feature-file output-filename-prefix \n";
	print "  Both input files are plain text and have NO header line.\n";
	print "  genomic-feature-file: the first four (4) columns are chrom, start, end, and category\n";
	exit(1);
}

sub uniq
{
	my %h;
	for(@_) {
		$h{$_} = 0;
	}
	return (sort keys %h);
}

