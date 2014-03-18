function [ lambda, wald_stat, p_val, snps ] = lin_reg( X, y, covar )
% Do a wald test to compute lin reg test statistic.  Variance parameter is
% estimated by REML.

[n, m] = size(X);
d = size(covar, 2) + 1;

wald_stat = zeros(1, m);

for i = 1:m
    % impute missing genotypes
    c = double(X(:, i));
    missing = c == 9;
    mu = mean(c(~missing));
    c(missing) = mu;
    
    Q = [c, covar];
    
    beta_hat = (Q'*Q)\(Q'*y);
    r = y - Q*beta_hat; % residuals
    sigma_e_sq = sum(r.^2)/(n - d);
    
    QQ_inv = inv(Q'*Q);
    wald_stat(i) = (beta_hat(1)^2/QQ_inv(1, 1))/sigma_e_sq;
end

p_val = 1 - chi2cdf(wald_stat, 1);
lambda = lambda_GC(wald_stat);
[~, snps] = sort(wald_stat, 'descend');

end
