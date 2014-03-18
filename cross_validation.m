function [ res ] = cross_validation(working, pheno_file, covar_file, ...
    top_k, n_folds, seed, use_PC, debug)
%{
Calculates CV log-likelihood over n_folds using the top_k SNPs.  Returns
the negative log-likelihood averaged over n_folds.
    
    Assumes that the working directory has:
    geno.mat
    pcs.mat

    Calculated from preproccess.m

    pheno_file: phenotypes in one columns
    covar_file: covariates in white space separated columns
    top_k: # of SNPs to use in GRM
    n_folds: number of folds to use in CV
    seed: random seed to choose folds
    use_PC: use the PCs as covariates  
%}
    
if nargin < 8
    debug = false;
end

%% Load in phenotypes
f = fopen(pheno_file, 'r');
C = textscan(f, '%s');
fclose(f);

% If case/control convert to binary
if strcmp(C{1}{1}, 'Case') || strcmp(C{1}{1}, 'Control')
    y = double(strcmp(C{1}, 'Case'));
else
    y = str2double(C{1});
end

n = length(y);

missing = isnan(y) | y == -9 | y == 9;
y = y(~missing);

fprintf('Read in phenotypes for %d individuals.\n', length(y));

%% Load in covariates
f = fopen(covar_file, 'r');
C = fscanf(f, '%f');
fclose(f);

n_covar = length(C)/n;
covar = reshape(C, n_covar, n)';

fprintf('Reading %d covariates.\n', n_covar);

% add all ones to covar
covar = [ones(size(covar, 1), 1) covar];

if use_PC
    load(sprintf('%s/pcs.mat', working)); % stored in U
    covar = [covar U];
    U = [];
end

covar = covar(~missing, :);

%% Run CV
load(sprintf('%s/geno.mat', working)); % stored in geno
geno = geno(~missing, :);

% for debugging, make smaller
if debug
    J = 1:1000;
    n = size(geno, 1);
    I = [1:50, n-50:n];
    res = cv_helper(geno(I, J), y(I), covar(I, :), top_k, n_folds, seed);
else
    res = cv_helper(geno, y, covar, top_k, n_folds, seed);
end

end