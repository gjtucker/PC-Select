function [ X, X_null, X_causal, y, I ] = sim_strat_and_causal_data( n, ...
    m, m_causal, m_null, p, f_st, trait_diff, sigma_e )
% If f_st = 0, then no population stratification

if f_st == 0 && trait_diff ~= 0
    fprintf('Possible error, F_st = 0, but trait_diff is nonzero');
end

m_tot = m + m_causal + m_null;

% Sample X
% sample randomly from uniform [0.1, 0.9]
x = rand(1, m_tot)*0.8 + 0.1;

if f_st == 0
    X = bsxfun(@lt, rand(n, m_tot), x) + bsxfun(@lt, rand(n, m_tot), x);
else
    % Sample population allele frequency using f_st equation
    x_1 = betarnd(x*(1 - f_st)/f_st, (1 - x)*(1 - f_st)/f_st);
    x_2 = betarnd(x*(1 - f_st)/f_st, (1 - x)*(1 - f_st)/f_st);

    % Sample realization of SNP matrix
    X = bsxfun(@lt, rand(n/2, m_tot), x_1) + bsxfun(@lt, rand(n/2, m_tot), x_1);
    X = [X; bsxfun(@lt, rand(n/2, m_tot), x_2) + bsxfun(@lt, rand(n/2, m_tot), x_2)];
end

X_causal = X(:, 1:m_causal);
X_null = X(:, (m_causal+1):(m_causal + m_null));
X(:, 1:(m_causal + m_null)) = [];

x = mean(X)/2;
W = bsxfun(@minus, X, 2*x);
W = bsxfun(@rdivide, W, sqrt(2*x.*(1 - x)));

x = mean(X_causal)/2;
W_causal = bsxfun(@minus, X_causal, 2*x);
W_causal = bsxfun(@rdivide, W_causal, sqrt(2*x.*(1 - x)));

% Sample causal markers
I = rand(m, 1) < p;
grm_effects = zeros(m, 1);

grm_effects(I) = randn(nnz(I), 1)*sqrt(0.5/nnz(I));
test_effects = randn(m_causal, 1)*sqrt(0.5/m_causal);

y = [W_causal, W] * [test_effects; grm_effects] + ...
    trait_diff * [ones(n/2, 1); zeros(n/2, 1)] + ...
    sigma_e * randn(n, 1);

end

