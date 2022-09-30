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
	push @{$snp{$c[0]}},[@c[3,1]];
}
close IN;

for (keys %snp)
{
	@{$snp{$_}} = sort {$a->[0] <=> $b->[0]} @{$snp{$_}};
}

print "There are ", scalar(@snp_arr), " SNPs in $snp_pos_file.\n";

# step 2
# get genomic feature
open IN,$gf_file or die "Cannot open file $gf_file: $!\n";
my %gf;
while(<IN>)
{
	chomp;
	my @c = split /\s+/;
	next unless (@c == 4);
	push @{$gf{$c[0]}},[@c[1..3]];
}
close IN;

for (keys %gf)
{
	@{$gf{$_}} = sort {$a->[0] <=> $b->[0]} @{$gf{$_}};
}

print "Completed reading $gf_file\n";

# step 3
# assign SNPs to genomic features with single pass
my %snp2gf;
my %useful_gf;
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
				push @{$snp2gf{$snp_id}}, $gf_c;
				$useful_gf{$gf_c} = 0;
			}
		}
	}
}
my @categories = sort keys %useful_gf;

if(scalar(keys %snp2gf) == 0) {
	print "Something went wrong: no SNPs are within genomic features.\n";
	print "No snp info file generated\n";
	exit(1);
}
print scalar(keys %snp2gf), " SNPs are assigned into ", scalar(@categories), " genomic features.\n";

# step 4
# generate output file
my %cat2idx;
for(0..$#categories) {
	$cat2idx{$categories[$_]} = $_;
}
my @null = ('') x scalar(@categories);
open OUT, ">$output_prefix.snp_info.csv";
print OUT join(",", "SNP", "intercept", @categories), "\n";
for (@snp_arr) {
	my @out = @null;
	if(defined $snp2gf{$_}) {
		my @cur_cats = &uniq(@{$snp2gf{$_}});
		for my $c (@cur_cats) {
			$out[ $cat2idx{$c} ] = 1;
		}
	}
	print OUT join(",", $_, 1, @out), "\n";
}

print "$output_prefix.snp_info.csv generated\n";

sub print_usage()
{
	print "Run the program:\n";
	print "  perl make_snp_info.pl plink-bim-file genomic-feature-file output-filename-prefix \n";
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

