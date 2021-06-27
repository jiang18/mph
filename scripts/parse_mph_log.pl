use strict;
use warnings;

@ARGV == 2 or die "2 arguments needed: log-filename and output-filename\n";
my ($log, $out) = @ARGV;

my @info;
my ($trt, $chol, $lu, $diff, $iter) = (0) x 5;
my @keys = qw(trt chol lu diff iter);
open IN,$log or die "Could not find $log: $!\n";
while(<IN>) {
  if(/trait_name with ARG (.*)/) {
    $trt = $1;
  }
  if(/completed normally/) {
    $chol ++;
  }
  if(/Switched to LU/) {
    $lu ++;
  }
  if(/diff for all VCs = ([\d\.]+)/) {
    $diff = $1;
  }
  if(/Completed MINQUE iteration (\d+)/) {
    $iter ++;
  }
  if(/Analysis finished/) {
    push @info, {trt=>$trt, chol=>$chol, lu=>$lu, diff=>$diff, iter=>$iter};
    $chol = 0;
    $lu = 0;
    $iter = 0;
  }
}
close IN;

open OUT,">$out" or die "Could not create $out: $!\n";
print OUT join(",",@keys), "\n";
for my $h (@info) {
  print OUT join(",", map{$h->{$_}} @keys),"\n";
}
close OUT;

