% Set random number generator 
rng default;

n_samples = 100;

%% Set up models
cutoffs = [10, 30, 100, 300, 1000, 3000, 10000];

models = cell(1, 1);

models{1}.name = 'Linear';
models{1}.m = Lin_model(false);

models{1}.name = 'Linear w/ PCs';
models{1}.m = Lin_model(true);

models{3}.name = 'LMM';
models{3}.m = MLM_model(false);

models{4}.name = 'FaST-LMM Select';
models{4}.m = Fast_lmm_autoselect_model(cutoffs, false);

models{5}.name = 'PC-Select';
models{5}.m = Fast_lmm_autoselect_model(cutoffs, true);

%% Set up data generators
m = 10000;
m_null = 500;
m_causal = 0;
p = 0;
sigma_e = 1;
f_st = 0.05;
trait_diff = 0.25;
n = 1000;

data_generators = cell(1, 1);

% No causal + pop strat
data_generators{1}.name = 'no causal + pop strat';
data_generators{1}.func = @(worker_id) sim_strat_and_causal_data( n, ...
    m, m_causal, m_null, p, f_st, trait_diff, sigma_e );

sigma_e = 0;

% p = 0.05 + pop strat
m_causal = 50;
p = 0.05;

data_generators{2}.name = 'p = 0.05 + pop strat';
data_generators{2}.func = @(worker_id) sim_strat_and_causal_data( n, ...
    m, m_causal, m_null, p, f_st, trait_diff, sigma_e );

% p = 0.005 + pop strat
m_causal = 50;
p = 0.005;

data_generators{3}.name = 'p = 0.005 + pop strat';
data_generators{3}.func = @(worker_id) sim_strat_and_causal_data( n, ...
    m, m_causal, m_null, p, f_st, trait_diff, sigma_e );

% p = 0.05
m_causal = 50;
p = 0.05;
f_st = 0;
trait_diff = 0;

data_generators{4}.name = 'p = 0.05';
data_generators{4}.func = @(worker_id) sim_strat_and_causal_data( n, ...
    m, m_causal, m_null, p, f_st, trait_diff, sigma_e );

% p = 0.005
m_causal = 50;
p = 0.005;
f_st = 0;
trait_diff = 0;

data_generators{5}.name = 'p = 0.005';
data_generators{5}.func = @(worker_id) sim_strat_and_causal_data( n, ...
    m, m_causal, m_null, p, f_st, trait_diff, sigma_e );
%}
    
%% Run the pipeline
res = struct();
res.models = models;
res.data_generators = data_generators;
res.runs = cell(length(data_generators), 1);

for i = 1:length(data_generators)
    res.runs{i} = run_pipeline(data_generators{i}.func, models, n_samples);
end

plot_figure_1(res);
