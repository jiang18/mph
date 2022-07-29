---
title: Prediction
layout: template
filename: prediction.md
---

## Empirical best linear unbiased predictions (EBLUPs)
```
mph --pred --mq_file milk.chr --output milk.chr
```
MPH computes EBLUPs using the output of --minque and outputs them to a file with a suffix of .mq.blup.csv. For genomic partitioning, EBLUPs are the estimates of direct genomic values. 
