function [ lambda, wald_stat, p_val ] = fast_lmm( X, X_test, y, covar, worker_id )

base_filename = sprintf('geno%d', worker_id);

write_to_t_file(sprintf('../working/sim_geno_sim_pheno/%s', base_filename), ...
    X, X_test, y, covar);

cmd = sprintf(['../FaSTLMM.207.Linux/Linux_MKL/fastlmmc ', ...
            '-out ../working/sim_geno_sim_pheno/%s.out.txt ', ...
            '-tfile ../working/sim_geno_sim_pheno/%s_test ', ...
            '-tfileSim ../working/sim_geno_sim_pheno/%s_cov ', ...
            '-pheno ../working/sim_geno_sim_pheno/%s_pheno.txt ', ...
            '-mpheno 1 ', ...
            '-covar ../working/sim_geno_sim_pheno/%s_covar.txt'], ...
            base_filename, base_filename, base_filename, base_filename, base_filename);

system(cmd);

% Process output
output = importdata(sprintf('../working/sim_geno_sim_pheno/%s.out.txt', base_filename));

% Reorder the output by the actual SNP ordering
[~, I] = sort(output.data(:, 1));

wald_stat = output.data(I, 13);
p_val = output.data(I, 5);

lambda = lambda_GC(wald_stat);

end

