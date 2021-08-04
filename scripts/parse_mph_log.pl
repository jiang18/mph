use strict;
use warnings;

@ARGV == 2 or die "2 arguments needed: log-filename and output-filename\n";
my ($log, $out) = @ARGV;

my @info;
my @keys = qw(trt rcond grad_norm iter diff newton reject);
my ($trt, $rcond, $grad_norm, $iter, $diff, $newton, $reject);
open IN,$log or die "Could not find $log: $!\n";
while(<IN>) {
  if(/trait_name with ARG (.*)/) {
    $trt = $1;
  }
  if(/number of H \= ([\d\-\.e]+)/) {
    $rcond = $1;
  }
  if(/Gradient norm \= ([\d\-\.e]+)/) {
    $grad_norm = $1;
  }
  if(/Started MINQUE iteration (\d+)/) {
    $iter ++;
  }
  if(/delta logLL \= ([\d\-\.e]+)/) {
    $diff = $1;
  }
  if(/descent \+ ([\d\-\.e]+)/) {
    $newton = $1;
  }
  if(/Reject the step/) {
    $reject ++;
  }
  
  if(/Analysis finished/) {
    push @info, {trt=>$trt, rcond=>$rcond, grad_norm=>$grad_norm, iter=>$iter, diff=>$diff, newton=>$newton, reject=>$reject};
    $iter = 0;
    $reject = 0;
  }
}
close IN;

open OUT,">$out" or die "Could not create $out: $!\n";
print OUT join(",",@keys), "\n";
for my $h (@info) {
  print OUT join(",", map{$h->{$_}} @keys),"\n";
}
close OUT;

