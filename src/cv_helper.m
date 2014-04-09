function [ total, rank_order ] = cv_helper(geno, y, covar, top_k, n_folds, seed, ...
    rank_order)

if nargin < 7
    recompute_rank_order = true;
    rank_order = zeros(n_folds, size(geno, 2)); % store for future use
else
    recompute_rank_order = false;
end

% Create folds
rng(seed, 'twister');
I_folds = randi(n_folds, length(y), 1);

% Do CV
total = 0;

for k = 1:n_folds
    fprintf('Starting fold %d\n', k);
    % Order SNPs by linear regression p-value
    I_train = I_folds ~= k;
    
    if recompute_rank_order
        [~, wald_stat] = lin_reg(geno(I_train, :), y(I_train), covar(I_train, :));
        [~, rank_order(k, :)] = sort(wald_stat, 'descend');
        wald_stat = [];
    end
    
    GRM = make_GRM_int8(geno(:, rank_order(k, 1:top_k)));
    [U, S] = svd(GRM(I_train, I_train));
    S = diag(S);
    fprintf('Finished SVD\n');
    
    [log_delta, sigma_g_sq] = opt_delta(U'*y(I_train), ...
        U'*covar(I_train, :), S);
    
    total = total + eval_nLL(I_train, y, ...
            GRM, covar, U, S, log_delta, sigma_g_sq);
        
    GRM = [];
    U = [];
    S = [];
end

end

function [log_delta, sigma_g_sq] = opt_delta(Uy, Ucovar, S)

% Optimize delta
grid = linspace(-5, 5, 100);
model = @(x) neg_log_lik(x, Ucovar, Uy, S);
coarse_log_delta = grid_search(model, grid);

% Fine optimization
log_delta = fminbnd(model, coarse_log_delta - 0.1, ...
        coarse_log_delta + 0.1);
    
[~, ~, ~, sigma_g_sq] = neg_log_lik(log_delta, Ucovar, Uy, S);
end

% negative log likelihood up to a constant and scaling
function nLL = eval_nLL(I_train, y, GRM, covars, U, S, log_delta, sigma_g_sq)

delta = exp(log_delta);
V_inv = 1./(S + delta);

Ucovars_train = U'*covars(I_train, :);
Uy_train = U'*y(I_train);

alpha = (Ucovars_train'*(bsxfun(@times, V_inv, Ucovars_train))) \ ...
    (Ucovars_train' * (V_inv .* Uy_train));

% Remove fixed effects
WWT = GRM(~I_train, I_train) * U;

y = y - covars*alpha;
y_pred = WWT*(V_inv.*Uy_train - V_inv.*(Ucovars_train*alpha));
y_diff = y(~I_train) - y_pred;

Sigma = sigma_g_sq*(GRM(~I_train, ~I_train) + delta*eye(nnz(~I_train)) - ...
    WWT * (bsxfun(@times, V_inv, WWT')));
R = chol(Sigma);

nLL = (log(2*pi)*nnz(~I_train) + 2*sum(log(diag(R))) + y_diff'*(R\(R'\y_diff)))/2;

end