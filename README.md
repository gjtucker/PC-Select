# PC-Select
Selects the number of SNPs to include in the GRM by CV log-likelihood.

To run, first convert files to eigenstrat format (http://www.hsph.harvard.edu/alkes-price/software/).

This example assumes that you are in the directory of .m files and you have a data directory with the genotypes in eigenstrat format.  The phenotype is data/phen.txt and covariates are data/covars.txt.

1) Run a preprocessor to generate working files

```matlab
n_pcs = 5;
preprocess('data/MS.geno', 'data', n_pcs);
```

Creates data/GRM.mat, data/pcs.mat.

2) Run CV code 

```matlab
top_k_choices = [10, 100, 1000, 10000];
res = zeros(length(top_k_choices), 1);

for i = 1:length(top_k_choices)
    res(i) = cross_validation('data', 'data/phen.txt', ...
        'data/covars.txt', top_k_choices(i), 10, 1, true);
end

[~, I] = min(res);
top_k_choices(I)
```
