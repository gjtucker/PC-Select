function W = normalize_genotypes(X)
x = mean(X)/2;
W = bsxfun(@minus, X, 2*x);
W = bsxfun(@rdivide, W, sqrt(2*x.*(1 - x)));
end