function preprocess(eigenstrat_genotypes, working, n_pcs)
% Read in .geno EIGENSTRAT format file and outputs genotype matrix in 
% mat format (GRM.mat) and PCs in mat format (pcs.mat) 
% 
% working: directory to output GRM and PCs
% n_pcs: number of pcs to save

if nargin < 2 
    eigenstrat_genotypes = '../data/WTCCC2/MS_ALL';
    working = '../working/WTCCC2';
end

if nargin < 3
    n_pcs = 5;
end

%% Load in genotypes
geno_file_name = sprintf('%s/geno.mat', working);
if exist(geno_file_name, 'file') == 2
    load(geno_file_name);
else
    geno = geno2mat(eigenstrat_genotypes, geno_file_name);
end

%% Compute GRM
GRM = make_GRM_int8(geno);
save(sprintf('%s/GRM.mat', working), 'GRM', '-v7.3');

%% Calc PCs
[U, ~] = svd(GRM);
GRM = []; % save memory
U = U(:, 1:n_pcs); % save top n_pcs

save(sprintf('%s/pcs.mat', working), 'U', '-v7.3');

end
