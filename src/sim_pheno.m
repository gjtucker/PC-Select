function [X, X_null, X_causal, y, I ] = sim_pheno(X, X_null, X_causal, ...
    ancestry_trait, p, multiplier, sigma_e )

W = normalize_genotypes(X);
W_causal = normalize_genotypes(X_causal);

[n, m] = size(X);
m_causal = size(X_causal, 2);

% Sample causal markers
I = rand(m, 1) < p;
grm_effects = zeros(m, 1);

grm_effects(I) = randn(nnz(I), 1)*sqrt(0.5/nnz(I));
test_effects = randn(m_causal, 1)*sqrt(0.5/m_causal);

y = [W_causal, W] * [test_effects; grm_effects] + ...
    ancestry_trait * multiplier + ...
    sigma_e * randn(n, 1);

end