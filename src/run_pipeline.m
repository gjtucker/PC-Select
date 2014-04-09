function res = run_pipeline(data_generator_func, models, n_samples)

res = cell(n_samples, 1);

parfor i = 1:n_samples
    res{i} = run_pipeline_helper(data_generator_func, models, i);
end
