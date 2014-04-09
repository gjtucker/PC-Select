function res = combine_results(working_path, n_samples, m_null, m_causal)

res = struct();
res.models = create_models();
res.data_generators = create_data_generators();

res.runs = cell(length(res.data_generators), 1);

for i = 1:length(res.data_generators)
    res.runs{i} = cell(n_samples, 1);
    for j = 1:n_samples
        res.runs{i}{j} = cell(length(res.models), 1);
        for k = 1:length(res.models)
            res.runs{i}{j}{k} = read_results( ...
                sprintf('%s/results/res_data_%d_sample_%d_model_%d', ...
                working_path, i, j, k), m_null, m_causal);
        end
    end
end

end

function res = read_results(filename, m_null, m_causal)

% Process output
output = importdata(sprintf('%s.out.txt', filename));

% Reorder the output by the actual SNP ordering
[~, I] = sort(output.data(:, 1));

wald_stat = output.data(I, 13);

res = struct();
res.wald_stat_null = wald_stat(1:m_null);
res.wald_stat_causal = wald_stat(m_null+1:end);

assert(m_null + m_causal == length(wald_stat));

end

% these are specific to the simulation that we ran!
function models = create_models()

models = cell(1, 1);

models{1}.name = 'Linear';
models{2}.name = 'Linear w/ PCs';
models{3}.name = 'LMM';
models{4}.name = 'FaST-LMM Select';
models{5}.name = 'PC-Select';

end

function data_generators = create_data_generators()

data_generators = cell(1);

data_generators{1}.name = 'p = 0.05 + pop strat';
data_generators{2}.name = 'p = 0.005 + pop strat';
data_generators{3}.name = 'p = 0.05';
data_generators{4}.name = 'p = 0.005';

end