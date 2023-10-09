use strict;
use warnings;

@ARGV == 1 or die "annot index needed\n";

my $annot = shift;

my $funct = "funct.snp_info.csv";
open IN,$funct;
$_=<IN>;
chomp;
my @head = split /,/;

$annot < @head or die "# of annotations: ",scalar(@head)-1,"\n";
my $annot_name = $head[$annot];

open OUT,">extract.txt";
while(<IN>) {
	chomp;
	my @c = split /,/;
	if(defined $c[$annot] and $c[$annot] ne '') {
		print OUT "$c[0]\n";
	}
}
close IN;
close OUT;

system("plink --bfile allchrs --extract extract.txt --make-bed --out subset --chr-set 30 --threads 10");

system("mph --verbose --compute_grm --binary_genotype subset --snp_info $funct --snp_weight $annot_name --output $annot_name --num_threads 10");
