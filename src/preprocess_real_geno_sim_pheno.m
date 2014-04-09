rng default;

data_prefix = '../data/MS_ALL'; % .bim
working = '../data'; % geno.mat
working_output = '../working/5k_individuals';

% Read in genotype matrix
geno_file_name = [working, '/geno.mat'];
load(geno_file_name); % stored in geno

% Read in chromosome information
f = fopen(sprintf('%s.bim', data_prefix), 'r');
C = textscan(f, '%d %s %d %d %s %s');
chrom = C{1};
C = [];
fclose(f);

% Subsample genotypes
[N, M] = size(geno);

n = 5000;
m = 50000;
m_causal = 250;
m_null = 500;

I = datasample(1:N, n, 'Replace', false);

% Sample test SNPs from chrom 1 and GRM SNPs from other chromosomes.  This
% is to avoid proximal contamination.
J = [datasample(find(chrom ~= 1 & chrom ~= 2), m, 'Replace', false); ...
       datasample(find(chrom == 1), m_causal, 'Replace', false); ...
       datasample(find(chrom == 2), m_null, 'Replace', false)];

geno = double(geno(I, J));

% Impute missing snps
for j = 1:size(geno, 2)
    c = geno(:, j);
    missing = c == 9;
    mu = mean(c(~missing));
    c(missing) = mu;
    geno(:, j) = c;
end

%% Save the subsampled genotypes
X = geno(:, 1:m);
X_causal = geno(:, m+1:m + m_causal);
X_null = geno(:, m+m_causal+1:end);

W = normalize_genotypes(X);
[U, S] = svd(W*W'/m);

% Extract PC1
ancestry_trait = U(:, 1);
ancestry_trait = ancestry_trait/std(ancestry_trait);

mkdir(sprintf('%s/genotypes', working_output));
mkdir(sprintf('%s/phenotypes', working_output));
mkdir(sprintf('%s/covariates', working_output));
mkdir(sprintf('%s/results', working_output));
mkdir(sprintf('%s/sh_src', working_output));
mkdir(sprintf('%s/logs', working_output));

% Write out data for Matlab to read in
save(sprintf('%s/genotypes/subsampled_genotypes.mat', working_output), ...
    'X', 'X_causal', 'X_null', 'ancestry_trait');

% Write out data for Fast LMM Select
base_filename_real_geno = 'geno';
write_to_t_file(sprintf('%s/genotypes/%s', working_output, base_filename_real_geno), ...
    X, [X_null, X_causal], ones(n, 1), U(:, 1:5));

write_covar(sprintf('%s/covariates/%s', working_output, 'pcs_covar.txt'), U(:, 1:5));
write_covar(sprintf('%s/covariates/%s', working_output, 'nocovar_covar.txt'), U(:, end+1:end));


