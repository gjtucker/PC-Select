# Accompanying code for PC-Select

To run analysis, first set up the directory:
* Run setup.sh 
* Download FaSTLMM Select from (http://research.microsoft.com/en-us/um/redmond/projects/mscompbio/fastlmm/) and put it in this directory.
* (Optional) The real genotypes and phenotypes are made available by the WTCCC2 (http://www.wtccc.org.uk/ccc2/).  To run those analyses, acquire the data from WTCCC2 and convert to eigenstrat format ((http://www.hsph.harvard.edu/alkes-price/software/).

All of the following commands should be run from the src directory.

## Simulated genotypes and phenotypes

```matlab
run_sim_geno_sim_pheno
```

Runs the simulation in Table 1.  Briefly, creates 100 simulated data sets with 1000 individuals and 10000 markers with p = 0.05, 0.005 markers causal and with or without population stratification.  See manuscript for details.

## Real genotypes and simulated phenotypes

Running this example requires real genotype data in plink binary format in the data directory.  First, convert the data into eigenstrat format. 

1) Run a preprocessor to generate working files

```matlab
n_pcs = 5;
preprocess('data/MS.geno', 'data', n_pcs);
preprocess_real_geno_sim_pheno
generate_data_and_runner('../working/5k_individuals', 1, 200)
```

2) Then run

```
sh ../working/5k_individuals/sh_src/run.sh
```

Runs the simulation in Table 2.  Briefly, creates 200 simulated data sets with 5000 individuals and 50000 markers (subsampled from real genotypes) with p = 0.05, 0.005 markers causal and with or without population stratification.  See manuscript for details.

3) To plot results

```matlab
res = combine_results('../working/5k_individuals/results', 200, 500, 250);
plot_figure_1(res);
```

## Cross-validation on large data sets

Selects the number of SNPs to include in the GRM by CV log-likelihood.

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
## Real genotypes and real phenotypes

For the WTCCC2 MS data, we found that cross-validation selected all markers for both PC-Select and FaST-LMM Select.  We used GCTA (http://www.complextraitgenomics.com/software/gcta/) to compute the association statistics in this case.  We ran the following to do that

```
gcta64 --bfile data/MS_ALL --autosome --make-grm --out data/MS_ALL --thread-num 8
gcta64 --grm data/MS_ALL --pca 5 --out data/MS_ALL
gcta64 --mlma-loco --bfile data/MS_ALL --pheno data/phen.txt --out data/MS_ALL_pc --thread-num 8 --mlma-no-adj-covar --qcovar data/MS_ALL.eigenvec
gcta64 --mlma-loco --bfile data/MS_ALL --pheno data/phen.txt --out MS_ALL_no_pc --thread-num 8 --mlma-no-adj-covar
```

where data/MS_ALL.{bim,bam,fam} contain the genotypes and phen.txt contains the phenotypes.



