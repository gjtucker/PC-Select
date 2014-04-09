function res = run_pipeline_helper(data_generator_func, models, seed)
% Set random seed
rng(seed);
res = cell(length(models), 1);

[X, X_null, X_causal, y]  = feval(data_generator_func, seed);

m_null = size(X_null, 2);

for j = 1:length(models)
	[wald_stat, p_val, extra] = models{j}.m.run(X, [X_null, X_causal], y, seed);

	wald_stat_null = wald_stat(1:m_null);
	wald_stat_causal = wald_stat(m_null+1:end);

	p_val_null = p_val(1:m_null);
	p_val_causal = p_val(m_null+1:end);

	res{j} = struct('wald_stat_null', wald_stat_null, ...
	    'wald_stat_causal', wald_stat_causal, ...
	    'p_val_null', p_val_null', ...
	    'p_val_causal', p_val_causal', ...
	    'extra', extra);
end

end

