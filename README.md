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

Runs the simulation in Table 1.  Briefly, creates 100 simulated data sets with 1000 individuals and 10000 markers with p = 0.05, 0.005 markers causal and with or without population stratification.  Results are plotted.

## Real genotypes and simulated phenotypes

Running this example requires real genotype data in plink binary format in the data directory.  First, convert the data into eigenstrat format. 

Run a preprocessor to generate working files

```matlab
n_pcs = 5;
preprocess('data/MS.geno', 'data', n_pcs);
preprocess_real_geno_sim_pheno
generate_data_and_runner('../working/5k_individuals', 1, 200)
```

Then run

```
sh ../working/5k_individuals/sh_src/run.sh
```

Runs the simulation in Table 2.  Briefly, creates 200 simulated data sets with 5000 individuals and 50000 markers (subsampled from real genotypes) with p = 0.05, 0.005 markers causal and with or without population stratification.

To plot results

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
