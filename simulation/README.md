```sh
mkdir familinx
cd familinx
perl clean_familinx_step1.pl
perl clean_familinx_step2.pl
cd ..
```
```sh
perl subset_ped.pl 1
cd 1
../markersim
../genosim
perl ../aipl2plink.pl 10k
plink --file 10k --make-bed --out 10k
```


