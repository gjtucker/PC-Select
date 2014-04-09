function generate_data_and_runner(working_path, seed, n_samples, debug)

if nargin < 4
    debug = false;
end

% Set random number generator 
rng(seed);

%% Set up data generators
load(sprintf('%s/genotypes/subsampled_genotypes.mat', working_path)); % loads X, X_null, X_causal, ancestry_trait

trait_multiplier = 0.25;
sigma_e = 0;

data_generators = cell(1, 1);

% p = 0.05 + pop strat
p = 0.05;

data_generators{1}.name = 'p = 0.05 + pop strat';
data_generators{1}.func = @(worker_id) sim_pheno(X, X_null, X_causal, ...
    ancestry_trait, p, ...
    trait_multiplier, sigma_e);


% p = 0.005 + pop strat
p = 0.005;

data_generators{2}.name = 'p = 0.005 + pop strat';
data_generators{2}.func = @(worker_id) sim_pheno(X, X_null, X_causal, ...
    ancestry_trait, p, ...
    trait_multiplier, sigma_e);

% p = 0.05
p = 0.05;
trait_multiplier = 0;

data_generators{3}.name = 'p = 0.05';
data_generators{3}.func = @(worker_id) sim_pheno(X, X_null, X_causal, ...
    ancestry_trait, p, ...
    trait_multiplier, sigma_e);

% p = 0.005
p = 0.005;

data_generators{4}.name = 'p = 0.005';
data_generators{4}.func = @(worker_id) sim_pheno(X, X_null, X_causal, ...
    ancestry_trait, p, ...
    trait_multiplier, sigma_e);
   
%% Write scripts to run FastLMM
cutoffs = [10, 30, 100, 300, 1000, 3000, 10000, 30000, 50000];

fastlmmc_path = '../FaSTLMM.207.Linux/Linux_MKL/fastlmmc';
script_file = fopen(sprintf('%s/sh_src/run.sh', working_path), 'w');
genotype_path = sprintf('%s/genotypes/geno', working_path);
test_snps = genotype_path;
grm_snps = genotype_path;
no_covar = sprintf('%s/covariates/nocovar', working_path);
pcs = sprintf('%s/covariates/pcs', working_path);

for i = 1:length(data_generators)
    for j = 1:n_samples
        [~, ~, ~, y] = data_generators{i}.func(j); 
                
        % write y to file
        write_y(sprintf('%s/phenotypes/data_%d_sample_%d_pheno.txt', working_path, ...
            i, j), y);
        
        pheno = sprintf('%s/phenotypes/data_%d_sample_%d', working_path, ...
                         i, j);
        
        models = cell(5, 1);
        models{1} = linear_cmd(fastlmmc_path, sprintf('%s/results/res_data_%d_sample_%d_model_1', working_path, i, j), ...
            test_snps, pheno, no_covar);
        
        models{2} = linear_cmd(fastlmmc_path, sprintf('%s/results/res_data_%d_sample_%d_model_2', working_path, i, j), ...
            test_snps, pheno, pcs);

        models{3} = fast_lmm_cmd(fastlmmc_path, sprintf('%s/results/res_data_%d_sample_%d_model_3', working_path, i, j), ...
            test_snps, grm_snps, pheno, no_covar);
        
        models{4} = fast_lmm_autoselect_cmd(fastlmmc_path, sprintf('%s/results/res_data_%d_sample_%d_model_4', working_path, i, j), ...
            test_snps, grm_snps, pheno, no_covar, cutoffs);
        
        models{5} = fast_lmm_autoselect_cmd(fastlmmc_path, sprintf('%s/results/res_data_%d_sample_%d_model_5', working_path, i, j), ...
            test_snps, grm_snps, pheno, pcs, cutoffs);
        
        if debug
            fprintf(script_file, '%s && %s && %s\n', models{1}, models{2}, models{3});
            fprintf(script_file, '%s\n', models{4});
            fprintf(script_file, '%s\n', models{5});
        else
            prefix = sprintf('bsub -q short -W 11:59 -R "rusage[mem=15000]" -o %s/logs/data_%d_sample_%d_models_%s.txt', ...
                working_path, i, j, '123');
            fprintf(script_file, '%s ''%s && %s && %s''\n', prefix, models{1}, models{2}, models{3});
            
            prefix = sprintf('bsub -q short -W 11:59 -R "rusage[mem=15000]" -o %s/logs/data_%d_sample_%d_models_%s.txt', ...
                working_path, i, j, '4');
            fprintf(script_file, '%s ''%s''\n', prefix, models{4});
            
            prefix = sprintf('bsub -q short -W 11:59 -R "rusage[mem=15000]" -o %s/logs/data_%d_sample_%d_models_%s.txt', ...
                working_path, i, j, '5');
            fprintf(script_file, '%s ''%s''\n', prefix, models{5});
        end  
    end
end

fclose(script_file);

end

function cmd = linear_cmd(fast_lmm, out, test_snps, pheno, covar)

cmd = sprintf(['%s ', ...
    '-out %s.out.txt ', ...
    '-tfile %s_test -linreg ', ...
    '-pheno %s_pheno.txt ', ...
    '-mpheno 1 ', ...
    '-covar %s_covar.txt'], ...
    fast_lmm, out, test_snps, pheno, covar);

end

function cmd = fast_lmm_cmd(fast_lmm, out, test_snps, grm_snps, pheno, covar)

cmd = sprintf(['%s ', ...
            '-out %s.out.txt ', ...
            '-tfile %s_test ', ...
            '-tfileSim %s_cov ', ...
            '-pheno %s_pheno.txt ', ...
            '-mpheno 1 ', ...
            '-covar %s_covar.txt'], ...
            fast_lmm, out, test_snps, grm_snps, pheno, covar);
  
end

function cmd = fast_lmm_autoselect_cmd(fast_lmm, out, test_snps, grm_snps, ...
    pheno, covar, cutoffs)

cutoff_strings = cell(length(cutoffs), 1);
for i = 1:length(cutoffs)
    cutoff_strings{i} = int2str(cutoffs(i));
end

autoSelectSearchValues = ['"', sprintf('%s ', cutoff_strings{1:end-1}), ...
    cutoff_strings{end}, '"'];

cmd = sprintf(['%s ', ...
            '-autoSelect %s ', ...
            '-autoSelectSearchValues %s ', ...
            '-randomSeed 1 ', ...
            '-autoSelectFolds 10 ', ...
            '-tfileSim %s_cov ', ...
            '-pheno %s_pheno.txt ', ...
            '-mpheno 1 ', ...
            '-covar %s_covar.txt'], ...
            fast_lmm, out, autoSelectSearchValues, grm_snps, pheno, covar);
        
cmd = sprintf(['%s && %s ', ...
            '-out %s.out.txt ', ...
            '-tfile %s_test ', ...
            '-tfileSim %s_cov ', ...
            '-extractSim %s.snps.txt ', ...
            '-pheno %s_pheno.txt ', ...
            '-mpheno 1 ', ...
            '-covar %s_covar.txt'], ...
            cmd, fast_lmm, out, test_snps, grm_snps, out, pheno, covar);
end
